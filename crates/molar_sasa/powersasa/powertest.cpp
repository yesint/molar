#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include <power_diagram.h>
#include <testincludes.h>


using namespace POWER_DIAGRAM;

int main()
{
#ifdef _FLOATTEST
	//Specifying Number Type and Vector Type
	typedef float number_type;
	typedef Eigen::Vector3f Coord; //fake vector4->3f might be even faster
#else
	typedef double number_type;
	typedef Eigen::Vector3d Coord;
#endif
        typedef std::vector<Coord,Eigen::aligned_allocator<Coord> > PDVector;
	PDVector coords; //Vector of coordinates
	std::vector<number_type> weights; //Vector of weights
	std::vector<int> bond_to; //best guess for connectivity


	coords.push_back(Coord(1.0,2.0,3.0));
	weights.push_back(0.56789);
	bond_to.push_back(0);
	coords.push_back(Coord(5.0,3.0,3.0));
	weights.push_back(1.12);
	bond_to.push_back(0);
	coords.push_back(Coord(5.1,3.2,3.5));
	weights.push_back(1.23);
	bond_to.push_back(1);
	coords.push_back(Coord(1.3,3.2,1.5));
	weights.push_back(1.34);
	bond_to.push_back(2);
	coords.push_back(Coord(0.0,3.4,1.7));
	weights.push_back(1.45);
	bond_to.push_back(3);
	coords.push_back(Coord(5.4,3.4,1.7));
	weights.push_back(1.56);
	bond_to.push_back(4);
	coords.push_back(Coord(3.4,2.1,1.4));
	weights.push_back(1.67);
	bond_to.push_back(5);

	std::cout << "Testing Power Diagram, 7 verts: " << std::endl;
	PowerDiagram<number_type,Coord,3> pd = PowerDiagram<number_type,Coord,3>::create(coords.size(),coords.begin(),weights.begin(),bond_to.begin())
		.with_radiiGiven(1).with_calculate(1).with_cells(1).with_Warnings(0);
	bond_to.clear();
	parse("lucy25K",coords,weights);
	bond_to.push_back(0);
	const double scale=1;
	coords[0]*=scale;
	weights[0]*=scale;
	for (unsigned int i = 1; i < coords.size(); ++i)
	{
		coords[i]*=scale;
		weights[i]*=scale;
		bond_to.push_back(i-1);
	}
	std::cout << "Testing Power Diagram, bunny_set: " << std::endl;
	PowerDiagram<number_type,Coord,3> pd2 = PowerDiagram<number_type,Coord,3>::create(coords.size(),coords.begin(),weights.begin(),bond_to.begin())
		.with_radiiGiven(1).with_calculate(1).with_cells(1).with_Warnings(0);
	//pd2.update_coords(coords);


	return 0;
}
