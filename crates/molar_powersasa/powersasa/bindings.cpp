#include "power_sasa.h"

/*
template<class ValueType>
class Selection_container_it_t {
public:
    using value_type = ValueType ;
    using difference_type = int ;
    using pointer = ValueType*;
    using reference = ValueType& ;
    using iterator_category = std::forward_iterator_tag;

    Selection_container_it_t(Selection* sel, int n) {parent = sel; pos = n;}

    Selection_container_it_t operator++(int junk) {auto tmp = *this; ++pos; return tmp;}
    Selection_container_it_t& operator++() { ++pos; return *this; }
    reference operator*() const { return parent->xyz(pos); }
    pointer operator->() { return parent->xyz_ptr(pos); }
    bool operator==(const Selection_container_it_t& rhs) { return pos == rhs.pos && parent == rhs.parent; }
    bool operator!=(const Selection_container_it_t& rhs) { return pos != rhs.pos || parent != rhs.parent; }

    operator Selection_container_it_t<const Eigen::Vector3f>() const { return *this; }
private:
    Selection* parent;
    int pos;
};

class Selection_coord_container {
public:
    Selection_coord_container(Selection& sel): parent(&sel){}

    using iterator = Selection_container_it_t<Eigen::Vector3f>;
    using const_iterator = Selection_container_it_t<const Eigen::Vector3f>;

    iterator begin(){ return iterator(parent,0); }
    const_iterator begin() const{ return const_iterator(parent,0); }
    const_iterator cbegin() const{ return const_iterator(parent,0); }
    iterator end(){ return iterator(parent,parent->size()); }
    const_iterator end() const { return const_iterator(parent,parent->size()); }
    const_iterator cend() const { return const_iterator(parent,parent->size()); }

    int size() const {return parent->size();}
private:
    Selection* parent;
};
*/

float powersasa(float probe_r, float* area_per_atom,
    float *total_volume, float* volume_per_atom)
{
    /*
    // First obtain the VDW radii of all atoms in selection and add probe radius to them
    vector<float> radii(size());
    for(int i=0; i<size(); ++i) radii[i] = vdw(i) + probe_r;

    bool do_v = total_volume ? true : false;
    bool do_a_per_atom = area_per_atom ? true : false;
    bool do_v_per_atom = volume_per_atom ? true : false;

    // Create aux container which allows not to copy selection coordinates locally
    Selection_coord_container cont(const_cast<Selection&>(*this));

    // Call POWERSASA
    POWERSASA::PowerSasa<float,Eigen::Vector3f> ps(cont, radii, 1, 0, do_v || do_v_per_atom, 0);
    ps.calc_sasa_all();

    float v,surf=0.0;

    if(do_v || do_v_per_atom){
        if(volume_per_atom) volume_per_atom->resize(size());
        for(int i = 0; i < size(); ++i){
            v = ps.getVol()[i];
            if(do_v_per_atom) (*volume_per_atom)[i] = v;
            if(do_v) *total_volume += v;
        }
    }

    if(area_per_atom) area_per_atom->resize(size());
    for(int i = 0; i < size(); ++i){        
        v = ps.getSasa()[i];
        //cout << i << " " << v << endl;
        if(do_a_per_atom) (*area_per_atom)[i] = v;
        surf += v;
    }

    return surf;
    */
}
