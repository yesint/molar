#ifndef __TESTINCLUDES_H_
#define __TESTINCLUDES_H_

#include <string>
#include <sstream>
#include <fstream>

template <class PDFloat,class PDCoord>
void parse (std::string const& filename, std::vector<PDCoord> & coord_output, std::vector<PDFloat>& weight_output)
{
	coord_output.clear();
	weight_output.clear();
	std::string temp;
	std::ifstream input(filename.c_str());
	while (getline(input,temp))
	{
		std::stringstream instream(temp);
		PDCoord tempcoord;
		PDFloat tempweight;
		instream >> tempcoord.x() >> tempcoord.y() >> tempcoord.z() >> tempweight ;
		coord_output.push_back(tempcoord);
		weight_output.push_back(tempweight);
	}
}


#endif // __TESTINCLUDES_H_
