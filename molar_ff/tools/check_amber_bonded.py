#!/usr/bin/env python3
"""Validate parameter completion against parmchk2, then export OpenMM references.

See assets/amber/README.md for the pinned source and reference-build correction.
Python dependencies: openmm==8.6.1, parmed==4.3.1, numpy.
"""
import argparse
import json, pathlib, subprocess, os, itertools, math
import tempfile
import parmed, openmm
from openmm import unit
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--models', required=True, type=pathlib.Path)
parser.add_argument('--parmchk2', required=True, type=pathlib.Path)
parser.add_argument('--amber-home', required=True, type=pathlib.Path)
parser.add_argument('--output', type=pathlib.Path)
args = parser.parse_args()
root = pathlib.Path(__file__).resolve().parents[1]
output = args.output or root / 'tests/data/amber_bonded_reference.json'
records = json.loads(args.models.read_text())
workspace = tempfile.TemporaryDirectory(prefix='molar-amber-reference-')
work = pathlib.Path(workspace.name)
errors = []
fixtures = []
methods = {}
for record in records:
    name, ff = (record['name'], record['ff'])
    ac = (root / 'tests/data/gaff_ref/ac' / f'{name}.ac').read_text().splitlines()
    out = []
    idx = 0
    for line in ac:
        if line.startswith('ATOM'):
            line = line.rsplit(None, 1)[0] + '  ' + record['types'][idx]
            idx += 1
        out.append(line)
    infile = work / f'{name}-{ff}.ac'
    infile.write_text('\n'.join(out) + '\n')
    extra = work / f'{name}-{ff}.frcmod'
    subprocess.run([str(args.parmchk2.resolve()), '-i', str(infile), '-f', 'ac', '-o', str(extra), '-s', str(2 if ff == 'gaff2' else 1)], env={**os.environ, 'AMBERCLASSICHOME': str(args.amber_home.resolve()), 'AMBERHOME': str(args.amber_home.resolve())}, check=True, capture_output=True)
    pars = parmed.amber.AmberParameterSet(str(root / 'assets/amber' / f'{ff}.dat'), str(extra))
    base = parmed.amber.AmberParameterSet(str(root / 'assets/amber' / f'{ff}.dat'))
    system = openmm.System()
    [system.addParticle(12.0) for p in record['positions']]
    bondforce = openmm.HarmonicBondForce()
    angleforce = openmm.HarmonicAngleForce()
    properforce = openmm.PeriodicTorsionForce()
    improperforce = openmm.PeriodicTorsionForce()
    for group, force in enumerate([bondforce, angleforce, properforce, improperforce]):
        force.setForceGroup(group)
        system.addForce(force)
    for kind in ['bonds', 'angles', 'propers', 'impropers']:
        for term in record[kind]:
            ids = term['atoms']
            types = tuple((record['types'][i] for i in ids))
            method = term['method']
            methods[method] = methods.get(method, 0) + 1
            if kind == 'bonds':
                p = pars.bond_types[types]
                k = p.k * 418.4
                r = p.req * 0.1
                if abs(k - term['k']) > 4.2 or abs(r - term['r']) > 6e-05:
                    errors.append((name, ff, kind, types, method, 'k/r', term['k'], k, term['r'], r, term['matched']))
                original = base.bond_types.get(tuple(term['matched']))
                if original is not None:
                    k = original.k * 418.4
                    r = original.req * 0.1
                bondforce.addBond(*ids, r, 2 * k)
            elif kind == 'angles':
                p = pars.angle_types[types]
                k = p.k * 4.184
                theta = math.radians(p.theteq)
                if abs(k - term['k']) > 0.05 or abs(theta - term['theta']) > 2e-05:
                    errors.append((name, ff, kind, types, method, 'k/theta', term['k'], k, term['theta'], theta, term['matched']))
                original = base.angle_types.get(tuple(term['matched']))
                if method == 'EmpiricalAngle':
                    a, b, c = term['matched']
                    r1 = base.bond_types[a, b].req
                    r2 = base.bond_types[b, c].req
                    theta = math.radians(0.5 * (base.angle_types[a, b, a].theteq + base.angle_types[c, b, c].theteq))
                    elements = {}
                    coeff = {}
                    for line in (root / 'assets/amber/PARMCHK.DAT').read_text().splitlines():
                        f = line.split()
                        if f and f[0] == 'PARM' and (len(f) >= 7):
                            elements[f[1]] = int(f[6])
                    for line in (root / 'assets/amber' / ('PARM_BLBA_GAFF2.DAT' if ff == 'gaff2' else 'PARM_BLBA_GAFF.DAT')).read_text().splitlines():
                        f = line.split()
                        if len(f) >= 6 and f[1] == 'BA':
                            coeff[int(f[3])] = (float(f[4]), float(f[5]))
                    k = 143.9 * coeff[elements[a]][1] * coeff[elements[b]][0] * coeff[elements[c]][1] * math.exp(-2 * ((r1 - r2) / (r1 + r2)) ** 2) / ((r1 + r2) * math.sqrt(theta)) * 4.184
                elif original is not None:
                    k = original.k * 4.184
                    theta = math.radians(original.theteq)
                angleforce.addAngle(*ids, theta, 2 * k)
            else:
                if kind == 'propers':
                    ps = pars.dihedral_types.get(types, pars.dihedral_types.get(('X', types[1], types[2], 'X')))
                    force = properforce
                else:
                    candidates = []
                    for key, p in (base if method in ['Exact', 'Wildcard'] else pars).improper_periodic_types.items():
                        for outer in itertools.permutations([types[0], types[1], types[3]]):
                            perm = (outer[0], outer[1], types[2], outer[2])
                            if all((a == b or a == 'X' for a, b in zip(key, perm))):
                                candidates.append((key.count('X'), key, p))
                    ps = [min(candidates, key=lambda v: v[0])[2]] if candidates else None
                    force = improperforce
                if ps is None:
                    errors.append((name, ff, kind, types, 'NO REFERENCE'))
                    continue
                expected = sorted([(p.phi_k * 4.184, math.radians(p.phase), p.per) for p in ps], key=lambda p: (p[2], p[1]))
                actual = sorted(term['terms'], key=lambda p: (p[2], p[1]))
                if len(expected) != len(actual) or any((abs(a - b) > 0.005 for p, q in zip(actual, expected) for a, b in zip(p, q))):
                    errors.append((name, ff, kind, types, method, actual, expected, term['matched']))
                for k, phase, n in expected:
                    force.addTorsion(*ids, int(n), phase, k)
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('Reference'))
    context.setPositions(record['positions'])
    state = context.getState(getEnergy=True, getForces=True)
    energies = [context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) for i in range(4)]
    fixtures.append({'name': name, 'ff': ff, 'positions': record['positions'], 'atomic_numbers': record['atomic_numbers'], 'formal_charges': record['formal_charges'], 'connectivity': record['connectivity'], 'energies': energies, 'forces': state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer).tolist()})
    del context, integrator
print('METHODS', methods)
print('MISMATCHES', len(errors))
for e in errors[:35]:
    print(e)
if errors:
    raise SystemExit('Parameter comparison failed; reference fixture was not changed')
metadata = {'openmm_version': openmm.__version__, 'parmed_version': parmed.__version__, 'amber_classic_commit': '76b00b2e9d7dc41ab1ce7c7aa2fe52b3608ec002', 'parmchk2_reference_correction': 'improperparmnum2 = improperparmnum (parameter count, not molecule improper count)', 'direct_impropers': 'All peripheral permutations; use base-table terms before parmchk2 estimates', 'resolution_methods': methods, 'records': fixtures}
output.write_text(json.dumps(metadata, indent=2) + '\n')
workspace.cleanup()
print(f'Wrote {len(fixtures)} OpenMM references to {output}')
