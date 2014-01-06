/* utils.cpp
	Utility functions for Audapter
	Part of Audapter
	
	Shanqing Cai, 2013
*/

#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>

#include "ost.h"

using namespace std;

/* Utility functions */
/* Read lines from a file */
list<string> readLinesFromFile(const string & ostFN) {
	list<string> lines;

	ifstream inFile(ostFN);

	if (!inFile)
		return lines;

	while (inFile.good()) {
		string t_str;

		getline(inFile, t_str);
		lines.push_back(t_str);
	}

	return lines;
}

/* Strip string: strip white space at the beginning and the end */
string trimString(const string & str) {
	string str_out;
	size_t pos_0 = str.find_first_not_of(" \t");
	size_t pos_1 = str.find_last_not_of(" \t");

	if ( (pos_0 != string::npos) && (pos_1 >= pos_0) ) {
		str_out = string(str, pos_0, pos_1 - pos_0 + 1);
	}
	
	return str_out;
}

/* Split string at white spaces; save results in a list */
list<string> splitStringToList(const string & str) {
	list<string> strs_out; 

	istringstream iss(str);

	while (iss.good()) {
		string t_str;
		iss >> t_str;
		strs_out.push_back(t_str);
	}

	return strs_out;
}

/* Split string by whitespace */
std::vector<std::string> splitStringToVector(const std::string &s) {
  std::vector<std::string> tokens;
  std::istringstream iss(s);

  std::copy(std::istream_iterator<std::string>(iss), 
	    std::istream_iterator<std::string>(),
	    std::back_inserter<std::vector <std::string> > (tokens));

  return tokens;
}

/* Remove comments:
	Input arguments: 
		items:	input items
		cc:		comment character
*/
vector<string> removeComments(const vector<string> items, const char cc) {
	vector<string> items_nc;

	for(vector<string>::const_iterator it = items.cbegin(); 
		it != items.cend(); ++it) {
		if ( ((*it).size() > 0) && ((*it)[0] == cc) )
			break;

		items_nc.push_back(*it);
	}

	return items_nc;
}
