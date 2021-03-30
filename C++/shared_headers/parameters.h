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

Vector3d vectorExtract(std::string line){
	Vector3d out;
	size_t begin = line.find_last_of(":")+2U;
	std::string myvector = line.substr(begin, line.size());
	size_t coma1 = myvector.find_first_of(",");
	size_t coma2 = myvector.find_last_of(",");
	out(0) = std::stof(myvector.substr(      0U, coma1));
	out(1) = std::stof(myvector.substr(coma1+1U, coma2));
	out(2) = std::stof(myvector.substr(coma2+1U, myvector.size()));
	return out;
}

class Parameters {
public:

	int Shape;
	std::string meshname;
	bool add_wrinkles;
	bool add_ramp;
	bool Abaqus_output;
	bool Dune_output;
	int nb_plies;

	int cz_id, resin_id;
	bool isResin, isCZ;

	// Wrinkle parameters
	std::string WID;
	double wrinkleSize;
	Vector3d wrinklePos;
	Vector3d wrinkleDamp;
	double wrinkleOri;

	// Ramp parameters
	double rampSize;

	// Path to results
	std::string path_to_abaqus_result;
	std::string AbaqusOdbName;

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

		if (line.find("Shape")!=std::string::npos)
			Shape = std::stoi(extract(line));

		if (line.find("Resin_betw_plies")!=std::string::npos)
			isResin = std::stoi(extract(line));

		if (line.find("cohezive_elements")!=std::string::npos)
			isCZ = std::stoi(extract(line));

		if (line.find("np(i)")!=std::string::npos)
			nb_plies = std::stoi(extract(line));

		if (line.find("name(s)")!=std::string::npos)
			meshname = extract(line)+".msh";

		if (line.find("wrinkle")!=std::string::npos)
			add_wrinkles = std::stoi(extract(line));

		if (line.find("Ramp(b)")!=std::string::npos)
			add_ramp = std::stoi(extract(line));

		if (line.find("Abaqus_output")!=std::string::npos)
			Abaqus_output = std::stoi(extract(line));

		if (line.find("Dune_output")!=std::string::npos)
			Dune_output = std::stoi(extract(line));

		if (line.find("Rsize(f)")!=std::string::npos)
			rampSize = std::stof(extract(line));

		if (line.find("WID(s)")!=std::string::npos)
			WID = extract(line);

		if (line.find("Wsize(f)")!=std::string::npos)
			wrinkleSize = std::stof(extract(line));

		if (line.find("Wori(f)")!=std::string::npos)
			wrinkleOri = std::stof(extract(line));

		if (line.find("Wpos(f)")!=std::string::npos)
			wrinklePos = vectorExtract(line);

		if (line.find("Wdamp(f)")!=std::string::npos)
			wrinkleDamp = vectorExtract(line);

		if (line.find("Path2result(s)")!=std::string::npos)
			path_to_abaqus_result = extract(line);

		if (line.find("AbaqusOdbName(s)")!=std::string::npos)
			AbaqusOdbName = extract(line);

		std::getline(input, line);
	}

	cz_id = -1;
	resin_id = -1;
	if (isResin && isCZ){
		cz_id = nb_plies+2;
		resin_id = nb_plies+1;
	} else if (isResin && !isCZ) {
		resin_id = nb_plies+1;
	} else if (!isResin && isCZ) {
		cz_id = nb_plies+1;
	}
}


#endif /* end of include guard: PARAMETERS_H */
