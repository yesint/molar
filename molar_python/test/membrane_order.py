from pymolar import *

#--------------------------------------
class MyTask(AnalysisTask):
    def register_args(self,parser):
        parser.add_argument('--params')
        

    def pre_process(self):
        toml_str=open(self.args.params).read()
        self.memb = Membrane(self.src,toml_str)

        z0 = 5.6
        self.memb.reset_valid_lipids()
        self.memb.reset_groups()
        
        upper = []
        lower = []
        for lip in self.memb.lipids:
            z = lip.head_marker[2]
            if lip.sel[0].resname == 'POGL' and z>4.6 and z<6.6:
            #if lip.sel[0].resname == 'POGL' and z>z0:
               lip.valid = False                
            else:
               lip.valid = True

            if z > z0:
                upper.append(lip.id)
            else:
                lower.append(lip.id)
        
        
        self.memb.add_lipids_to_group("upper", upper);
        self.memb.add_lipids_to_group("lower", lower);


    def process_frame(self):
        self.memb.set_state(self.state)
        self.memb.compute()
        #print(self.memb.lipids[0].sel[0].pos)


    def post_process(self):
        self.memb.finalize()


if __name__ == "__main__":
    MyTask()
