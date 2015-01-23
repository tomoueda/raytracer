#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <stdio.h> 
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <CImg.h>

DECLARE_string(input);

void parseInputs();
void parseInputLine(std::vector<std::string> line);
std::vector<std::string> &split(const std::string &s, char delim, 
	std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);