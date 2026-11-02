#!/usr/bin/env python3
"""Generate independent OpenMM LJ/shifted-force Coulomb references.

Requires openmm==8.6.1, parmed==4.3.1, numpy, and the bonded test's
MOLAR_AMBER_EXPORT JSON. No Rust energy, force, or exclusion data are imported.
"""
import argparse
import json
import math
from pathlib import Path
import openmm
import parmed
from openmm import unit

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--models', required=True, type=Path)
parser.add_argument('--output', type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
output = args.output or root / 'tests/data/amber_nonbonded_reference.json'
inputs = []
for record in json.loads(args.models.read_text()):
    charges = [0.3 * math.sin(0.71 * i) for i in range(len(record['types']))]
    mean = sum(charges) / len(charges)
    inputs.append(dict(name=record['name'], ff=record['ff'], types=record['types'],
                       bonds=[b[:2] for b in record['connectivity']],
                       charges=[q - mean for q in charges], positions=record['positions']))
for ff in ['gaff', 'gaff2']:
    for distance in [0.25, 0.99, 1.0, 1.05, 1.19, 1.2, 1.21]:
        inputs.append(dict(name=f'pair_{distance}', ff=ff, types=['c3', 'o'],
                           bonds=[], charges=[0.4, -0.4], positions=[[0.,0.,0.], [distance,0.,0.]]))
    inputs.append(dict(name='one_four_beyond_cutoff', ff=ff, types=['c3']*4,
                       bonds=[[0,1],[1,2],[2,3]], charges=[0.2,-0.2,0.3,-0.3],
                       positions=[[0.,0.,0.],[0.5,0.1,0.],[1.,0.,0.],[1.5,0.,0.1]]))
    inputs.append(dict(name='zero_lj_hydrogen', ff=ff, types=['ho', 'o'], bonds=[],
                       charges=[0.2,-0.2], positions=[[0.,0.,0.],[0.5,0.,0.]]))

for record in inputs:
    parameters = parmed.amber.AmberParameterSet(str(root / 'assets/amber' / (record['ff'] + '.dat')))
    system = openmm.System()
    lj = openmm.NonbondedForce()
    lj.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)
    lj.setCutoffDistance(1.2)
    lj.setUseSwitchingFunction(True)
    lj.setSwitchingDistance(1.0)
    lj.setUseDispersionCorrection(False)
    coulomb = openmm.CustomNonbondedForce('k*q1*q2*(1/r-2/rc+r/(rc*rc))')
    coulomb.addGlobalParameter('k', 138.935457644382)
    coulomb.addGlobalParameter('rc', 1.2)
    coulomb.addPerParticleParameter('q')
    coulomb.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
    coulomb.setCutoffDistance(1.2)
    exceptions = openmm.CustomBondForce('k*qq/r')
    exceptions.addGlobalParameter('k', 138.935457644382)
    exceptions.addPerBondParameter('qq')
    for t, q in zip(record['types'], record['charges']):
        p = parameters.atom_types[t]
        system.addParticle(12.0)
        lj.addParticle(q, 2*p.rmin*0.1/(2**(1/6)), p.epsilon*4.184)
        coulomb.addParticle([q])
    # OpenMM independently generates shortest-path exclusions and 1-4 pairs.
    lj.createExceptionsFromBonds(record['bonds'], 1/1.2, 0.5)
    for i in range(lj.getNumExceptions()):
        a, b, qq, sigma, epsilon = lj.getExceptionParameters(i)
        coulomb.addExclusion(a, b)
        value = qq.value_in_unit(unit.elementary_charge**2)
        if value != 0:
            exceptions.addBond(a, b, [value])
        lj.setExceptionParameters(i, a, b, 0., sigma, epsilon)
    for i in range(lj.getNumParticles()):
        _, sigma, epsilon = lj.getParticleParameters(i)
        lj.setParticleParameters(i, 0., sigma, epsilon)
    lj.setForceGroup(0)
    coulomb.setForceGroup(1)
    exceptions.setForceGroup(1)
    for force in [lj, coulomb, exceptions]:
        system.addForce(force)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    context.setPositions(record['positions'])
    record['energies'] = [context.getState(getEnergy=True, groups={g}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) for g in [0,1]]
    record['forces'] = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer).tolist()
    del context, integrator
output.write_text(json.dumps(dict(openmm_version=openmm.__version__, parmed_version=parmed.__version__, cutoff_nm=1.2, lj_switch_nm=1.0, electrostatics='shifted force; unmodified scaled 1-4 exceptions', records=inputs), indent=2)+'\n')
print(f'Wrote {len(inputs)} reference cases to {output}')
