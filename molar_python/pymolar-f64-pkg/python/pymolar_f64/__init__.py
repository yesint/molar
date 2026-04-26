"""pymolar-f64: double-precision build of MolAR Python bindings.

This package mirrors :mod:`pymolar` but is compiled with the ``f64`` feature
of the underlying ``molar`` Rust library. All coordinates, masses, charges,
energies and time fields are 64-bit floats throughout the call stack and the
NumPy arrays exposed by the bindings have ``dtype=numpy.float64``.

The two wheels (`pymolar` and `pymolar_f64`) are co-installable in the same
environment because they live under distinct top-level packages and ship
distinct native modules (``pymolar.molar`` vs ``pymolar_f64.molar``).
"""

from .molar import *
import argparse
import logging

PBC_FULL = [True, True, True]
PBC_NONE = [False, False, False]
PBC_XY = [True, True, False]


def _process_suffix(s):
    if s == '':
        return (None, None)
    fr = None
    t = None
    if s[-2:].isnumeric():
        fr = int(s)
    elif s[-2:] == 'ps':
        t = int(s[:-2])
    elif s[-2:] == 'ns':
        t = int(s[:-2]) * 1000
    elif s[-2:] == 'us':
        t = int(s[:-2]) * 1000_000
    return (fr, t)


class AnalysisTask:
    """Base class for trajectory processing tasks (f64 build).

    See :class:`pymolar.AnalysisTask` for the full docstring; behavior is
    identical except that all numerical values are double precision.
    """

    def __init__(self):
        greeting()
        FORMAT = '[%(levelname)s] (%(name)s) %(message)s'
        logging.basicConfig(format=FORMAT)
        logging.getLogger().setLevel(logging.INFO)

        logging.info(f'Executing task "{type(self).__name__}"...')

        parser = argparse.ArgumentParser('molar_python trajectory processor (f64)')
        parser.add_argument('-f', '--files', nargs='+')
        parser.add_argument('--log', default=100, type=int)
        parser.add_argument('-b', '--begin', default='')
        parser.add_argument('-e', '--end', default='')
        parser.add_argument('--skip', default=1, type=int)
        parser.add_argument('--add-time', action="store_true")

        self.register_args(parser)
        self.args = parser.parse_args()

        if len(self.args.files) < 2:
            raise Exception('At least one trajectory file is required')

        self.top = FileHandler(self.args.files[0], 'r').read_topology()

        bfr, bt = _process_suffix(self.args.begin)
        efr, et = _process_suffix(self.args.end)

        self.consumed_frames = 0
        valid_frames = 0
        added_time = 0.0

        for trj_ind, trj_file in enumerate(self.args.files[1:]):
            logging.info(f'Processing trajectory "{trj_file}"...')
            self.trj_ind = trj_ind
            trj_handler = FileHandler(trj_file, 'r')

            if bfr:
                trj_handler.skip_to_frame(bfr)
            elif bt:
                trj_handler.skip_to_time(bt)

            for st in trj_handler:
                if efr and self.consumed_frames >= efr:
                    break
                if et and st.time + added_time > et:
                    break

                valid_frames += 1
                if (valid_frames - 1) % self.args.skip > 0:
                    continue

                st.time += added_time
                self.state = st

                if self.consumed_frames == 0:
                    self.src = System(self.top, self.state)
                    self.pre_process()
                else:
                    self.src.replace_state_deep(self.state)

                if self.consumed_frames % self.args.log == 0:
                    self.__log_time()

                self.consumed_frames += 1
                self.process_frame()

            if self.args.add_time:
                added_time += self.state.time

        self.post_process()

    def __log_time(self):
        if self.state.time < 1000.0:
            t = f"{self.state.time} ps"
        elif self.state.time < 1000_000.0:
            t = f"{self.state.time / 1000.0} ns"
        else:
            t = f"{self.state.time / 1000_000.0} us"
        logging.info(f'At frame {self.consumed_frames}, time {t}')

    def register_args(self, parser):
        pass

    def pre_process(self):
        pass

    def process_frame(self):
        pass

    def post_process(self):
        pass
