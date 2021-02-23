#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <sstream>
#include <regex>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Geometry>

using std::string;
using namespace Eigen;


std::vector<int> other(int n) {
	std::vector<int> tmp;
	tmp.resize(6);
	
	//~ std::vector<int> vec = { 7, 3, 5, 1, 9 };
	
	if      (n==0) {
		tmp[0] = 1;
		tmp[1] = 2;
		tmp[2] = 3;
		tmp[3] = 5;
		tmp[4] = 6;
		tmp[5] = 7;
	} else if (n==1) {
		tmp[0] = 0;
		tmp[1] = 2;
		tmp[2] = 3;
		tmp[3] = 4;
		tmp[4] = 6;
		tmp[5] = 7;
	} else if (n==2) {
		tmp[0] = 0;
		tmp[1] = 1;
		tmp[2] = 3;
		tmp[3] = 4;
		tmp[4] = 5;
		tmp[5] = 7;
	} else if (n==3) {
		tmp[0] = 0;
		tmp[1] = 1;
		tmp[2] = 2;
		tmp[3] = 4;
		tmp[4] = 5;
		tmp[5] = 6;
	}
	
	return tmp;
}

int diag(int n){
	int val;
	
	if (n==0) {val= 2;}
	else if (n==1) {val= 3;}
	else if (n==2) {val= 0;}
	else if (n==3) {val= 1;}
	
	return val;
}

int right(int n){
	int val;
	
	if (n==0) {val= 1;}
	else if (n==1) {val= 2;}
	else if (n==2) {val= 3;}
	else if (n==3) {val= 0;}
	
	return val;
}

int left(int n){
	int val;
	
	if (n==0) {val= 3;}
	else if (n==1) {val= 0;}
	else if (n==2) {val= 1;}
	else if (n==3) {val= 2;}
	
	return val;
}

std::vector<int> splitngetINT(std::string& line) {
	std::vector<int> out;
	out.resize(0);
	auto start = 0U;
	auto end = line.find(" ");
	while (end != std::string::npos) {
		out.push_back(std::stoi(line.substr(start, end - start)));
		start = end + 1U;
		end = line.find(" ", start);
	}
	out.push_back(std::stoi(line.substr(start, end)));
	return out;
}

bool isnul(VectorXf& v) {
	bool nul=true;
	float epsi = 0.00001;
	for(int i=0; i<v.size();i++) {if(v(i)>epsi) {nul=false;break;}}
	return nul;
}

bool issame(Vector3f old, Vector3f nouveau) {
	bool same=true;
	float epsi = 0.001;
	for(int i=0; i<3;i++) {
		if(std::abs(old(i)-nouveau(i))>epsi) {
			same=false;
			break;
		}
	}
	return same;
}

bool issame2(Vector3f old, Vector3f nouveau) {
	bool same=false;
	float epsi = 0.001;
	if(std::abs(old(0)-nouveau(0))<epsi && std::abs(old(1)-nouveau(1))<epsi && std::abs(old(2)-nouveau(2))<epsi) {
		same=true;
	}
	return same;
}

void strip_basename(const string& filename, std::string& pre, std::string& suf) {
	std::smatch m;
	std::regex e("(_)");

	std::regex_search(filename,m,e);
	pre = m.prefix().str();
	suf = m.suffix().str();
}

float ReverseFloat( const float inFloat ) {
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

unsigned short Reverse16bit( const unsigned short inFloat ) {
   unsigned short retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[1];
   returnFloat[1] = floatToConvert[0];

   return retVal;
}

std::string int2string(int f) {
  std::ostringstream os;
  os << f;
  return os.str();
}

std::string float2string(float f) {
  std::ostringstream os;
  os << f;
  return os.str();
}

#endif /* end of include guard: UTILS_H */
