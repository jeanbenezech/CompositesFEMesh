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
	out(0) = std::stod(myvector.substr(      0U, coma1));
	out(1) = std::stod(myvector.substr(coma1+1U, coma2));
	out(2) = std::stod(myvector.substr(coma2+1U, myvector.size()));
	return out;
}

Vector2d ExtractV2d(std::string line){
	Vector2d out;
	size_t begin = line.find_last_of(":")+2U;
	std::string myvector = line.substr(begin, line.size());
	size_t coma = myvector.find_first_of(",");
	out(0) = std::stod(myvector.substr(      0U, coma));
	out(1) = std::stod(myvector.substr(coma+1U, myvector.size()));
	return out;
}

class Parameters {
public:

	int Shape;
	std::string meshname;
	int add_wrinkles;
	bool add_ramp;
	bool Abaqus_output;
	bool Dune_output;
	int nb_plies;

	int cz_id, resin_id;
	bool isResin, isCZ;

	std::vector<Vector2d> StackSeq;

	// Wrinkle parameters
	std::string WID;
	std::vector<double> wrinkleSize;
	std::vector<Vector3d> wrinklePos;
	std::vector<Vector3d> wrinkleDamp;
	std::vector<double> wrinkleOri;

	// Ramp parameters
	double rampSize;
	double R;
	Vector3d StartEndinZdir;

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
		std::cout << "Error: Cannot open file" << filename << std::endl;
	} else {
		std::cout << "Reading " << filename << std::endl;
	}

	std::string line {};
	std::getline(input, line);

	StackSeq.resize(0);
	wrinkleSize.resize(0);
	wrinklePos.resize(0);
	wrinkleDamp.resize(0);
	wrinkleOri.resize(0);

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

		if (line.find("R(f)")!=std::string::npos)
			R = std::stod(extract(line));

		if (line.find("StartEndinZdir(f)")!=std::string::npos)
			StartEndinZdir = vectorExtract(line);

		if (line.find("Abaqus_output")!=std::string::npos)
			Abaqus_output = std::stoi(extract(line));

		if (line.find("Dune_output")!=std::string::npos)
			Dune_output = std::stoi(extract(line));

		if (line.find("Rsize(f)")!=std::string::npos)
			rampSize = std::stod(extract(line));

		if (line.find("WID(s)")!=std::string::npos)
			WID = extract(line);

		if (line.find("Wsize")!=std::string::npos){
			double tmp = std::stof(extract(line));
			wrinkleSize.push_back(tmp);
		}

		if (line.find("Wori")!=std::string::npos){
			double tmp = std::stof(extract(line));
			wrinkleOri.push_back(tmp);
		}

		if (line.find("Wpos")!=std::string::npos){
			Vector3d tmp = vectorExtract(line);
			wrinklePos.push_back(tmp);
		}

		if (line.find("Wdamp")!=std::string::npos){
			Vector3d tmp = vectorExtract(line);
			wrinkleDamp.push_back(tmp);
		}

		if (line.find("Path2result(s)")!=std::string::npos)
			path_to_abaqus_result = extract(line);

		if (line.find("AbaqusOdbName(s)")!=std::string::npos)
			AbaqusOdbName = extract(line);

		if (line.find("(f,f)")!=std::string::npos){
			Vector2d tmp = ExtractV2d(line);
			StackSeq.push_back(tmp);
		}

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
