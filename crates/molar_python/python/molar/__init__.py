from .molar import *
import argparse
import logging

PBC_FULL = [True,True,True]
PBC_NONE = [False,False,False]
PBC_XY = [True,True,False]

def _process_suffix(s):
    if s=='':
        return (None,None)
    fr = None
    t = None
    if s[-2:].isnumeric():
        #No suffix, frames provided
        fr = int(s)
    elif s[-2:] == 'ps':
        t = int(s[:-2])
    elif s[-2:] == 'ns':
        t = int(s[:-2])*1000
    elif s[-2:] == 'us':
        t = int(s[:-2])*1000_000
    return (fr,t)


class AnalysisTask:
    def __init__(self):
        greeting()
        FORMAT = '[%(levelname)s] (%(name)s) %(message)s'
        logging.basicConfig(format=FORMAT)
        logging.getLogger().setLevel(logging.INFO)

        logging.info(f'Executing task "{type(self).__name__}"...')

        parser = argparse.ArgumentParser('molar_python trajectory processor')
        parser.add_argument('-f','--files',nargs='+')
        parser.add_argument('--log',default=100,type=int)
        parser.add_argument('-b','--begin',default='')
        parser.add_argument('-e','--end',default='')
        parser.add_argument('--skip',default=1,type=int)

        # Register user-supplied arguments
        self.register_args(parser)
        # Parse arguments
        self.args = parser.parse_args()
        
        if len(self.args.files) < 2:
            raise Exception('At least one trajectory file is required')
        
        self.top = FileHandler(self.args.files[0]).read_topology()

        bfr,bt = _process_suffix(self.args.begin)
        efr,et = _process_suffix(self.args.end)

        # Read trajectories and call process_frame on each frame
        self.consumed_frames = 0
        begin_skipped = False
        valid_frames = 0
        
        for trj_file in self.args.files[1:]:
            logging.info(f'Processing trajectory "{trj_file}"...')
            trj_handler = _ParTrajReader(trj_file)
            
            if not begin_skipped:
                if bfr:
                    trj_handler.skip_to_frame(bfr)
                elif bt:
                    trj_handler.skip_to_time(bt)
                begin_skipped = True
            
            while True:
                # Read next frame
                st = trj_handler.next_state()
                # If None returned exit
                if st == None:
                    break

                # See if end is reached
                if efr and self.consumed_frames >= efr:
                    break
                if et and st.time > et:
                    break
                
                # We have a valid frame, 
                valid_frames += 1
                
                # see if we need to skip a frame
                if (valid_frames-1) % self.args.skip > 0:
                    continue

                self.state = st
                
                if self.consumed_frames == 0:
                    self.src = Source(self.top,self.state)
                    # Call pre-processing
                    self.pre_process()
                else:
                    self.src.set_state(self.state)
                        
                if self.consumed_frames % self.args.log == 0:
                    self.__log_time()

                self.consumed_frames += 1
                # User supplied process
                self.process_frame()

        # Call post-process
        self.post_process()


    def __log_time(self):
        if self.state.time < 1000.0:
            t = f"{self.state.time} ps"
        elif self.state.time < 1000_000.0:
            t = f"{self.state.time/1000.0} ns"
        else:
            t = f"{self.state.time/1000_000.0} us"
        logging.info(f'At frame {self.consumed_frames}, time {t}')


    def register_args(self,parser):
        pass


    def pre_process(self):
        pass


    def process_frame(self):
        pass


    def post_process(self):
        pass
