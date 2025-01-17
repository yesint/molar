from .molar import *
import argparse

class AnalyisTask:
    def __init__(self):
        parser = argparse.ArgumentParser('Trajectory processor')
        parser.add_argument('-f','--files',nargs='+')
        parser.add_argument('--log',default=100,type=int)
        parser.add_argument('-b','--begin',default=0,type=int)
        parser.add_argument('-e','--end',default=-1,type=int)
        args = parser.parse_args()
        for a in args.files:
            print(a)

        if len(args.files) < 2:
            raise Exception('At least one trajectory file is required')
        
        self.top = FileHandler(args.files[0]).read_topology()

        # Read trajectories and call process_frame on each frame
        self.consumed_frames = 0
        for trj_file in args.files[1:]:
            print(f'Processing trajectory "{trj_file}"...')
            trj_handler = FileHandler(trj_file)
            for st in trj_handler:
                self.state = st
                if self.consumed_frames == 0:
                    self.src = Source(self.top,self.state)
                    # Call pre-processing
                    self.pre_process()
                else:
                    self.src.set_state(self.state)
                        
                if self.consumed_frames % args.log == 0:
                    self.log_time()

                self.consumed_frames += 1
                self.process_frame()

        # Call post-process
        self.post_process()


    def log_time(self):
        if self.state.time < 1000.0:
            t = f"{self.state.time} ps"
        elif self.state.time < 1000_000.0:
            t = f"{self.state.time/1000.0} ns"
        else:
            t = f"{self.state.time/1000_000.0} us"
        print(f'Frame {self.consumed_frames}, time {t}')


    def pre_process(self):
        pass


    def process_frame(self):
        pass


    def post_process(self):
        pass
