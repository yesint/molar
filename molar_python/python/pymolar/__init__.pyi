"""Package API; coordinate arrays use this package's compiled precision."""
from .molar import *
import argparse

PBC_FULL: list[bool]
PBC_NONE: list[bool]
PBC_XY: list[bool]

class AnalysisTask:
    """Construction parses CLI arguments and runs the complete trajectory task."""
    args: argparse.Namespace
    top: Topology
    state: State
    src: System
    consumed_frames: int
    trj_ind: int
    def __init__(self) -> None: ...
    def register_args(self, parser: argparse.ArgumentParser) -> None: ...
    def pre_process(self) -> None: ...
    def process_frame(self) -> None: ...
    def post_process(self) -> None: ...
