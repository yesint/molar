#include "pteros/core/mol_file.h"
#include "pteros/core/format_recognition.h"
#include "pteros/core/pteros_error.h"
//---------------------------------------------------

using namespace std;
using namespace pteros;

int main(int argc, char** argv){
    try{
        //vector<molfile_plugin_t*> plugins;
        System sys;
        Frame fr;
        sys.frame_append(fr);

        boost::shared_ptr<Mol_file> h = io_factory("sel.pdb",'r');
        Mol_file_content c;
        c.structure = true;
        c.coordinates = true;
        h->read(&sys,&sys.Frame_data(0),c);

        /*
        boost::shared_ptr<Mol_file> h2 = io_factory("sel.pdb",'r');
        sys.frame_append(fr);
        h2->read(&sys,&sys.Frame_data(1),c);
        */

        Selection sel(sys,"all");
        cout << sel.size() << endl;
        cout << sel.Name(0) << endl;        
        cout << sys.Box(0) << endl;
        cout << sel.XYZ(0).transpose() << endl;
        //sel.set_frame(1);
        cout << sel.XYZ(0).transpose() << endl;              

        h.reset();

        // Test writing
        h = io_factory("res.pdb",'w');
        sel.set_frame(0);
        h->write(sel,c);
        h.reset();

        // Test DCD write
        h = io_factory("test.dcd",'w');
        c.structure = false;
        c.coordinates = false;
        c.trajectory = true;
        for(int i=0; i<sys.num_frames(); ++i){
            sel.set_frame(i);
            h->write(sel,c);
        }
        h.reset();

        // Test DCD read
        h = io_factory("test.dcd",'r');
        while(true){
            sys.frame_append(fr);
            bool ok = h->read(&sys,&sys.Frame_data(sys.num_frames()-1),c);
            cout << "frame" << endl;
            if(!ok){
                sys.frame_delete(sys.num_frames()-1);
                break;
            }
        }
        h.reset();

        // Test TRR write
        h = io_factory("test.xtc",'w');
        cout << sys.num_frames() << endl;
        for(int i=0; i<sys.num_frames(); ++i){
            cout << "Frame " << i << endl;
            sel.set_frame(i);
            h->write(sel,c);
        }
        h.reset();

        // Test TRR read
        h = io_factory("test.xtc",'r');
        int k =0;
        while(true){
            sys.frame_append(fr);
            bool ok = h->read(&sys,&sys.Frame_data(sys.num_frames()-1),c);
            cout << "frame " << k << endl;
            if(!ok){
                cout << "end on frame " << k << endl;
                sys.frame_delete(sys.num_frames()-1);
                break;
            }
            ++k;
        }
        h.reset();

        cout << "Testing load" << endl;
        // Test new System.load
        System s1("sel.pdb");
        s1.load("res.pdb");
        s1.load("test.dcd");
        s1.load("test.trr");
        s1.load("test.xtc");

    } catch(Pteros_error e) {
        cout << e.what() << endl;
    }
}
