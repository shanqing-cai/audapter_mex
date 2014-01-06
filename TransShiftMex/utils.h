/* utils.h
	Utility functions for Audapter
	Part of Audapter
	
	Shanqing Cai, 2013
*/


#ifndef UTILS_H
#define UTILS_H

/* String and I/O utilities */
/* Read lines from a file */
std::list<std::string> readLinesFromFile(const std::string & ostFN);

/* Strip string: strip white space at the beginning and the end */
std::string trimString(const std::string & str);

/* Split string at white spaces; save results in a list */
std::list<std::string> splitStringToList(const std::string & str);

/* Split string by whitespace */
std::vector<std::string> splitStringToVector(const std::string &s);

#endif