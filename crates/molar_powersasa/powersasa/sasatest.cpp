#include <stdio.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include <power_sasa.h>
#include <testincludes.h>

// g++ -o sasatest -I. sasatest.cpp

int main(int argc, char* argv[])
{
        if (argc != 2)
        {
                std::cout << "Usage: " << argv[0] << " <input file>" << std::endl;
                std::exit(1);
        }
   
  	//Specifying Number Type and Vector Type
	//typedef float number_type;
	//typedef Eigen::Vector3f Coord; //fake vector4->3f might be even faster
//#ifdef _FLOATTEST
	typedef float Scalar;
	typedef Eigen::Vector3f Coord;
	//typedef std::array<float,3> Coord;
	std::cout << "Precision: float" << std::endl;
//#else
	//typedef double Scalar;
	//typedef Eigen::Vector3d Coord;
	//std::cout << "Precision: double" << std::endl;
//#endif

	std::vector<Coord> coords; //Vector of coordinates
	std::vector<Scalar> weights; //Vector of weights
/*
	coords.push_back(Coord(1.0,2.0,3.0));
	weights.push_back(1.0);
	coords.push_back(Coord(5.0,3.0,3.0));
	weights.push_back(1.0);
	coords.push_back(Coord(5.1,3.2,3.5));
	weights.push_back(1.0);
	coords.push_back(Coord(1.3,3.2,1.5));
	weights.push_back(1.0);
	coords.push_back(Coord(0.0,3.4,1.7));
	weights.push_back(1.0);
	coords.push_back(Coord(5.4,3.4,1.7));
	weights.push_back(1.4);
	coords.push_back(Coord(3.4,2.1,1.4));
	weights.push_back(1.4);

	std::cout << "Testing SASA, 7 Verts" << std::endl;
	SASA::PowerSasa<number_type,Coord> *ps;
	ps = new SASA::PowerSasa<number_type,Coord>(coords, weights, 0, 0, 0, 0);
	//ps->calc_sasa_all();
	delete ps;
*/
        parse(argv[1],coords,weights);
	//std::vector<Scalar>                 Sasa;
	//std::vector< std::vector<Coord> >   DSasa;
	//std::vector<Scalar>                 Vol;
	//std::vector<Coord>                  DVol;
	std::vector<Coord>                  force;
	
	//Sasa.resize(coords.size());
	//DSasa.resize(coords.size());
	//Vol.resize(coords.size());
	//DVol.resize(coords.size());
	force.resize(coords.size());
	      
	POWERSASA::PowerSasa<Scalar,Coord> *ps = 
		new POWERSASA::PowerSasa<Scalar,Coord>(coords, weights, 1, 1, 1, 1);
	ps->calc_sasa_all();
	
    /*
	for (unsigned int i = 0; i < coords.size(); ++i) {
	  force[i] = -(ps->getDSasa())[i][ps->NumOfNeighbours(i)];
	}  
	for (unsigned int i = 0; i < coords.size(); ++i)
	{
		for (unsigned int j = 0; j < ps->NumOfNeighbours(i); ++j)
		{
			force[ps->AtomNo(i, j)] -=  ps->getDSasa()[i][j];
			// force[ps.AtomNo(i, j)] -=  sigma[i] * ps->getDSasa()[i][j];
		}
	}
    */
	
    Scalar volume = 0.0, surf = 0.0;
	for (unsigned int i = 0; i < coords.size(); ++i)
	{
		printf("%4d sasa=%7.3lf vol=%7.3lf\n", i, (ps->getSasa())[i], (ps->getVol())[i]);
		//printf("%4d sasa=%7.3lf fx=%8.4lf fy=%8.4lf fz=%8.4lf v=%8.4lf dvx=%8.4lf dvy=%8.4lf dvz=%8.4lf\n",
		//	i, ps->getSasa()[i], force[i][0], force[i][1], force[i][2], 
		//	ps->getVol()[i], ps->getDVol()[i][0], ps->getDVol()[i][1], ps->getDVol()[i][2]);
		volume += (ps->getVol())[i];
        surf += (ps->getSasa())[i];
	}
	printf("volume=%lf\n", volume);
    printf("sasa  =%lf\n", surf);

        //for (unsigned int i = 0; i < coords.size(); ++i)
        //{
        //        std::cout << coords[i][0] << " " << coords[i][1] << " " << coords[i][2]
        //        << " " << weights[i] << " " << Sasa[i] << std::endl;
        //}	
	
        delete ps;

	return 0;
}
