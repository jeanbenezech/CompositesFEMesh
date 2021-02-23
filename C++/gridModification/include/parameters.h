#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "utils.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

std::string extract(std::string line) {
	std::string out;
	size_t index = line.find_last_of(":");
	out = line.substr(index+2U, line.size());
	return out;
}

class Parameters {
public:

	std::string meshname;
	bool add_wrinkles;
	bool Abaqus_output;
	bool Dune_output;
	int nb_plies;
	int cz_id;

	void read(const std::string& filename);

	//~Parameters();
// private:
};

// ~~~~~~~~~~~~~~ READ ~~~~~~~~~~~~~~

void Parameters::read(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
		std::cout << "Error: Cannot open file" << std::endl;
	} else {
		std::cout << "Reading " << filename << std::endl;
	}

	std::string line {};
	std::getline(input, line);

	while(!line.empty()){

		if (line.find("Resin_betw_plies")!=std::string::npos)
			cz_id = std::stoi(extract(line));

		if (line.find("np(i)")!=std::string::npos)
			nb_plies = std::stoi(extract(line));

		if (line.find("name(s)")!=std::string::npos)
			meshname = extract(line)+".msh";

		if (line.find("wrinkle")!=std::string::npos)
			add_wrinkles = std::stoi(extract(line));

		if (line.find("Abaqus_output")!=std::string::npos)
			Abaqus_output = std::stoi(extract(line));

		if (line.find("Dune_output")!=std::string::npos)
			Dune_output = std::stoi(extract(line));

		std::getline(input, line);
	}
}


#endif /* end of include guard: PARAMETERS_H */
