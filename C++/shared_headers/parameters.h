#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "utils.h"
#include "../mINI/src/mini/ini.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#define PI 3.14159265

std::string extract(std::string line) {
	std::string out;
	size_t index = line.find_last_of(":");
	out = line.substr(index+2U, line.size());
	return out;
}

Vector3d vectorExtract(std::string line, bool isIni=false){
	Vector3d out;
	size_t begin = line.find_last_of(":")+2U;
	if (isIni == true)
		begin = 0U;
	std::string myvector = line.substr(begin, line.size());
	size_t coma1 = myvector.find_first_of(",");
	size_t coma2 = myvector.find_last_of(",");
	out(0) = std::stod(myvector.substr(      0U, coma1));
	out(1) = std::stod(myvector.substr(coma1+1U, coma2));
	out(2) = std::stod(myvector.substr(coma2+1U, myvector.size()));
	return out;
}

Vector2d ExtractV2d(std::string line, bool isIni=false){
	Vector2d out;
	size_t begin = line.find_last_of(":")+2U;
	if (isIni == true)
		begin = 0U;
	std::string myvector = line.substr(begin, line.size());
	size_t coma = myvector.find_first_of(",");
	out(0) = std::stod(myvector.substr(      0U, coma));
	out(1) = std::stod(myvector.substr(coma+1U, myvector.size()));
	return out;
}

std::vector<double> Extractvector(std::string line, int size, bool isIni=false){
	std::vector<double> out(size);
	size_t begin = line.find_last_of(":")+2U;
	if (isIni == true)
		begin = 0U;
	std::string myvector_init = line.substr(begin, line.size());
	std::string myvector = line.substr(begin, line.size());
	auto init_string = 0U;
	for (int i=0; i<size; i++){
		size_t coma = myvector.find_first_of(",");
		auto extracted = myvector.substr(0U, coma);
		out[i] = std::stod(extracted);
		// std::cout << "out["<<i<<"]= " << out[i] << std::endl;
		myvector = myvector_init.substr(coma + 1U, myvector_init.size());
		myvector_init = myvector;
	}
	return out;
}

std::vector<int> ExtractvectorINT(std::string line, int size, bool isIni=false){
	std::vector<int> out(size);
	size_t begin = line.find_last_of(":")+1U;
	if (isIni == true)
		begin = 0U;
	std::string myvector_init = line.substr(begin, line.size());
	std::string myvector = line.substr(begin, line.size());
	auto init_string = 0U;
	for (int i=0; i<size; i++){
		// std::cout << "line= " << myvector << std::endl;
		size_t coma = myvector.find_first_of(",");
		// std::cout << "coma= " << coma << std::endl;
		auto extracted = myvector.substr(0U, coma);
		// std::cout << "extracted= " << extracted << std::endl;
		out[i] = std::stoi(extracted);
		// std::cout << "out["<<i<<"]= " << out[i] << std::endl;
		myvector = myvector_init.substr(coma + 1U, myvector_init.size());
		myvector_init = myvector;
	}
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
	std::vector<int> cz_ids;
	bool isResin, isCZ, isShell, writeTransformedMSH;

	bool GaussianThickness, CornerThickness;
	double sigma, length, ThicknessVar;

	// std::vector<Vector2d> StackSeq;
	std::vector<Vector3d> StackSeq;

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

	// RVE rotation parameters
	double angleRVE = 0.0;
	double AngleRotateFlangeR = 0.0;
	double AngleRotateFlangeL = 0.0;
	double interior_radius_RVE = 0.0;
	bool rotateRVE = false;
	bool RotateFlanges = false;
	std::string rotateAxis = "X";
	Vector2d rotateStartStop;
	Vector2d centerRot;

	// Geometric parameters
	double X, Y, Z, Height;

	// Path to results
	std::string path_to_abaqus_result;
	std::string AbaqusOdbName;

	bool do_test_delam = false;
	bool do_flatten = false;
	int SizeCentersZ, SizeCentersX;
	std::vector<double> centersZ;
	std::vector<double> centersX;
	

	// Rigid boundary element
	int RigidBoundary=0;
	int SizeIntervaldZ=2;
	std::vector<double> intervaldZ;

	void read(const std::string& filename);

	//~Parameters();
// private:
};

// ~~~~~~~~~~~~~~ READ ~~~~~~~~~~~~~~

void Parameters::read(const std::string& filename) {
	StackSeq.resize(0);
	wrinkleSize.resize(0);
	wrinklePos.resize(0);
	wrinkleDamp.resize(0);
	wrinkleOri.resize(0);

	// INI reading
	mINI::INIFile file(filename+".ini");
	mINI::INIStructure config;
	file.read(config);

	// GENERAL parameters
	meshname = config["GENERAL"].get("name(s)")+".msh";
	if (config["GENERAL"].has("Shell(b)")==true)
		isShell = std::stoi(config["GENERAL"].get("Shell(b)"));
	if (config["GENERAL"].has("Shape(i)")==true)
		Shape = std::stoi(config["GENERAL"].get("Shape(i)"));
	if (config["GENERAL"].has("Auto_Resin_betw_plies(b)")==true)
		isResin = std::stoi(config["GENERAL"].get("Auto_Resin_betw_plies(b)"));
	if (config["GENERAL"].has("cohezive_elements(b)")==true)
		isCZ = std::stoi(config["GENERAL"].get("cohezive_elements(b)"));
	if (config["GENERAL"].has("Abaqus_output(b)")==true)
		Abaqus_output = std::stoi(config["GENERAL"].get("Abaqus_output(b)"));
	if (config["GENERAL"].has("Path2result(s)")==true)
		path_to_abaqus_result = config["GENERAL"].get("Path2result(s)");
	if (config["GENERAL"].has("AbaqusOdbName(s)")==true)
		AbaqusOdbName = config["GENERAL"].get("AbaqusOdbName(s)");
	if (config["GENERAL"].has("writeTransformedMSH(b)")==true)
		writeTransformedMSH = std::stoi(config["GENERAL"].get("writeTransformedMSH(b)"));
	if (config["GENERAL"].has("Dune_output(b)")==true)
		Dune_output = std::stoi(config["GENERAL"].get("Dune_output(b)"));




	// TESTDELAM
	if (config["TESTDELAM"].has("test_delam(b)")==true)
		do_test_delam = std::stoi(config["TESTDELAM"].get("test_delam(b)"));
	if (config["TESTDELAM"].has("flatten(b)")==true)
		do_flatten = std::stoi(config["TESTDELAM"].get("flatten(b)"));
	if (config["TESTDELAM"].has("sizeCentersZ(i)")==true)
		SizeCentersZ = std::stoi(config["TESTDELAM"].get("sizeCentersZ(i)"));
	if (config["TESTDELAM"].has("sizeCentersX(i)")==true)
		SizeCentersX = std::stoi(config["TESTDELAM"].get("sizeCentersX(i)"));
	if (do_test_delam==true){
		if (config["TESTDELAM"].has("centersZ(f)")==true){
			centersZ = Extractvector(config["TESTDELAM"].get("centersZ(f)"), SizeCentersZ, true);
		}
		if (config["TESTDELAM"].has("centersX(f)")==true){
			centersX = Extractvector(config["TESTDELAM"].get("centersX(f)"), SizeCentersX, true);
		}
	}



	// General Geometry
	if (config["GEOMETRY"].has("np(i)")==true)
		nb_plies = std::stoi(config["GEOMETRY"].get("np(i)"));
	if (config["GEOMETRY"].has("X(f)")==true)
		X = std::stod(config["GEOMETRY"].get("X(f)"));
	if (config["GEOMETRY"].has("Y(f)")==true)
		Y = std::stod(config["GEOMETRY"].get("Y(f)"));
	if (config["GEOMETRY"].has("Z(f)")==true)
		Z = std::stod(config["GEOMETRY"].get("Z(f)"));
	if (config["GEOMETRY"].has("Height(f)")==true)
		Height = std::stod(config["GEOMETRY"].get("Height(f)"));
	if (config["GEOMETRY"].has("R(f)")==true)
		R = std::stod(config["GEOMETRY"].get("R(f)"));

	// StaqSeq
	StackSeq.resize(nb_plies);
	for(int iter_ply = 0; iter_ply<nb_plies; iter_ply++){
		if (config["STACKSEQ"].has("p"+std::to_string(iter_ply)+"(f,f,b)")==true)
			StackSeq[iter_ply] = vectorExtract(config["STACKSEQ"].get("p"+std::to_string(iter_ply)+"(f,f,b)"), true); // true : from config file, false: from txt older parameters file format
	}

	/// Grid transformation
	if (config["GEO-TRANSFORMATION"].has("Ramp(b)")==true)
		add_ramp = std::stoi(config["GEO-TRANSFORMATION"].get("Ramp(b)"));
	if (config["GEO-TRANSFORMATION"].has("StartEndinZdir(f)")==true)
		StartEndinZdir = vectorExtract(config["GEO-TRANSFORMATION"].get("StartEndinZdir(f)"), true);
	if (config["GEO-TRANSFORMATION"].has("Rsize(f)")==true)
		rampSize = std::stod(config["GEO-TRANSFORMATION"].get("Rsize(f)"));
	if (config["GEO-TRANSFORMATION"].has("RotateRVE(b)")==true)
		rotateRVE = std::stoi(config["GEO-TRANSFORMATION"].get("RotateRVE(b)"));
	if (config["GEO-TRANSFORMATION"].has("AngleRotateRVE(f)")==true)
		angleRVE = std::stod(config["GEO-TRANSFORMATION"].get("AngleRotateRVE(f)"));
	// if (config["GEO-TRANSFORMATION"].has("InteriorRadiusRVE(f)")==true)
	// 	interior_radius_RVE = std::stod(config["GEO-TRANSFORMATION"].get("InteriorRadiusRVE(f)"));
	rotateAxis = config["GEO-TRANSFORMATION"].get("RotateAxis(s)");

	bool need_default_rotate_start_stop=true;
	if (rotateAxis=="X"){
		rotateStartStop = {0., Z };
	} else if (rotateAxis=="Z"){
		rotateStartStop = {0., X };
	}
	if (config["GEO-TRANSFORMATION"].has("Rotate_start_end(f)")==true)
		rotateStartStop = ExtractV2d(config["GEO-TRANSFORMATION"].get("Rotate_start_end(f)"), true);
	if (config["GEO-TRANSFORMATION"].has("RotateFlanges(b)")==true)
		RotateFlanges = std::stoi(config["GEO-TRANSFORMATION"].get("RotateFlanges(b)"));
	if (config["GEO-TRANSFORMATION"].has("AngleRotateFlangeRi(f)")==true)
		AngleRotateFlangeR = std::stod(config["GEO-TRANSFORMATION"].get("AngleRotateFlangeRi(f)"));
	if (config["GEO-TRANSFORMATION"].has("AngleRotateFlangeLe(f)")==true)
		AngleRotateFlangeL = std::stod(config["GEO-TRANSFORMATION"].get("AngleRotateFlangeLe(f)"));

	if (config["GEO-TRANSFORMATION"].has("RigidBoundary(b)")==true)
		RigidBoundary = std::stoi(config["GEO-TRANSFORMATION"].get("RigidBoundary(b)"));
	if (config["GEO-TRANSFORMATION"].has("dZIntervalsize(i)")==true)
		SizeIntervaldZ = std::stoi(config["GEO-TRANSFORMATION"].get("dZIntervalsize(i)"));
	if (RigidBoundary==true){
		if (config["GEO-TRANSFORMATION"].has("dZInterval(f)")==true){
			intervaldZ = Extractvector(config["GEO-TRANSFORMATION"].get("dZInterval(f)"), SizeIntervaldZ, true);
			// for (auto value : intervaldZ)
			// 	std::cout << value << std::endl;
		}
	}
	// GAUSSIAN parameters
	if (config["GENERAL"].has("GaussianThickness(b)")==true)
		GaussianThickness = std::stoi(config["GENERAL"].get("GaussianThickness(b)"));
	if (config["GAUSSIAN"].has("Sigma(d)")==true)
		sigma = std::stod(config["GAUSSIAN"].get("Sigma(d)"));
	if (config["GAUSSIAN"].has("Length(d)")==true)
		length = std::stod(config["GAUSSIAN"].get("Length(d)"));

	// std::cout << "length: " << length << std::endl;

	// CORNER THICKNESS parameters
	if (config["GENERAL"].has("CornerThickness(b)")==true)
		CornerThickness = std::stoi(config["GENERAL"].get("CornerThickness(b)"));
	if (config["GEO-TRANSFORMATION"].has("ThicknessVar(d)")==true)
		ThicknessVar = std::stod(config["GEO-TRANSFORMATION"].get("ThicknessVar(d)"));
	// std::cout << "CornerThickness: " << CornerThickness << std::endl;
	// std::cout << "ThicknessVar: " << ThicknessVar << std::endl;
	// WRINKLES
	if (config["WRINKLES"].has("wrinkle(i)")==true)
		add_wrinkles = std::stoi(config["WRINKLES"].get("wrinkle(i)"));
	WID = config["WRINKLES"].get("WID(s)");

	wrinkleSize.resize(add_wrinkles);
	wrinkleOri.resize(add_wrinkles);
	wrinklePos.resize(add_wrinkles);
	wrinkleDamp.resize(add_wrinkles);
	for(int iter_w = 0; iter_w<add_wrinkles; iter_w++){
		if (config["WRINKLES"].has("Wsize"+std::to_string(iter_w)+"(f)")==true)
			wrinkleSize[iter_w] = std::stod(config["WRINKLES"].get("Wsize"+std::to_string(iter_w)+"(f)"));
		if (config["WRINKLES"].has("Wori"+std::to_string(iter_w)+"(f)")==true)
			wrinkleOri[iter_w] = std::stod(config["WRINKLES"].get("Wori"+std::to_string(iter_w)+"(f)"));
		if (config["WRINKLES"].has("Wpos"+std::to_string(iter_w)+"(f)")==true)
			wrinklePos[iter_w] = vectorExtract(config["WRINKLES"].get("Wpos"+std::to_string(iter_w)+"(f)"), true);
		if (config["WRINKLES"].has("Wdamp"+std::to_string(iter_w)+"(f)")==true)
			wrinkleDamp[iter_w] = vectorExtract(config["WRINKLES"].get("Wdamp"+std::to_string(iter_w)+"(f)"), true);
	}

	/// Ajustment on parameters choices to avoid incompatibility
	if (rotateRVE){
		if (angleRVE>0){
			interior_radius_RVE = (rotateStartStop(1)-rotateStartStop(0)) / (angleRVE*PI / 180.);
		} else {
			std::cout << "the angle of rotate RVE must be positive" << std::endl; 
			exit( 0 );
		}
	}


	// Abaqus CZ ids
	int cz_ids_size;
	if (config["GENERAL"].has("Abaqus_cohezive_ids_size(i)")==true){
		cz_ids_size = std::stoi(config["GENERAL"].get("Abaqus_cohezive_ids_size(i)"));
		// std::cout << "cz_ids_size: " << cz_ids_size << std::endl;
		cz_ids.resize(cz_ids_size);
		cz_ids= ExtractvectorINT(config["GENERAL"].get("Abaqus_cohezive_ids(i)"), cz_ids_size);
		// for (auto value : cz_ids)
		// 	std::cout << value << std::endl;
	}



	cz_id = -1;
	resin_id = -1;
	if (isResin==1 && isCZ==1){
		cz_id = nb_plies+2;
		resin_id = nb_plies+1;
	} else if (isResin==1 && !isCZ==1) {
		resin_id = nb_plies+1;
	} else if (!isResin==1 && isCZ==1) {
		cz_id = nb_plies+1;
	}

	if (Shape==0 || Shape==2) { // Cspar: force this parameter to be false
		rotateRVE= false;
	}

}


#endif /* end of include guard: PARAMETERS_H */
