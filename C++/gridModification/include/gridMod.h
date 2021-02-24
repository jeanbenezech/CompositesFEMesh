#pragma once

#include "mesh.h"
#include "element.h"
#include "parameters.h"

#include <cmath>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Geometry>

using namespace Eigen;

class GridTransformation{
	public:

	// Global gemoetry
	float z1, z2, z3, z4;
	float delta_max;
	float zmin, zmax, ymin, ymax, xmin, xmax;
	float ymid;

	// Wrinkles parameters
	Vector3f defect_location, damping, center_rot;
	float wrinkleOri;
	float ref_lenght;
	float defect_size;

	// General parameters
	int verbosity = 0;

	void initialise(Mesh& m, Parameters& param);
	void ramp(Vector3f& point);
	void setWrinklesParameter(Parameters& param);
	void Y_alligned_wrinkles(Vector3f& point);
	void wrinkles(Vector3f& point);

	private:
	// float epsilon = 0.1;
};

void GridTransformation::initialise(Mesh& m, Parameters& param) {
	zmax=0;
	zmin=m.vertices(2, 0);
	ymin=0;
	ymax=m.vertices(1, 0);
	xmin=0;
	xmax=m.vertices(0, 0);

	for(int node=0; node < m.Nb_vertices(); ++node) {
		if (m.vertices(2, node)>zmax){ zmax = m.vertices(2, node);}
		if (m.vertices(2, node)<zmin){ zmin = m.vertices(2, node);}
		if (m.vertices(1, node)>ymax){ ymax = m.vertices(1, node);}
		if (m.vertices(1, node)<ymin){ ymin = m.vertices(1, node);}
		if (m.vertices(0, node)>xmax){ xmax = m.vertices(0, node);}
		if (m.vertices(0, node)<xmin){ xmin = m.vertices(0, node);}
	}

	if (verbosity>0){
		std::cout << "Xmin : " << xmin << std::endl;
		std::cout << "Xmax : " << xmax << std::endl;
		std::cout << "Ymin : " << ymin << std::endl;
		std::cout << "Ymax : " << ymax << std::endl;
		std::cout << "Zmin : " << zmin << std::endl;
		std::cout << "Zmax : " << zmax << std::endl;
	}

	float deltaz = (zmax-zmin) / 16.0;
	z1 = zmin + 2.0*deltaz;
	z2 = zmin + 7.0*deltaz;
	z3 = zmin + 9.0*deltaz;
	z4 = zmin + 14.0*deltaz;

	delta_max = param.rampSize;

	ymid = (ymax-ymin)/2.0;

	// std::cout << "Ymid : " << ymid << std::endl;
}

void GridTransformation::ramp(Vector3f& point){
	float ymid = (ymax+ymin)/2.0;

	if (point(2)>z1 && point(2)<z4){ // section 1, 2, 3

		if (point(2)>=z2 && point(2)<=z3){ // section 2

			float initial_0 = point(0);
			float decrease_0 = initial_0 * delta_max / xmax; // or xmid?
			point(0) -= decrease_0;

			float initial_1 = point(1);
			if(initial_1 > ymid) {
				float decrease_1 = (delta_max/(ymax-ymid))* (initial_1 - ymid);
				point(1) -= decrease_1;
			} else {
				float increase_1 = (delta_max/(ymid-ymin))* (ymid - initial_1);
				point(1) += increase_1;
			}

		} else if (point(2)<z2){ // section 1

			float local_delta = (delta_max / (z2-z1)) * point(2) + ((z1*delta_max)/(z1-z2));

			float initial_0 = point(0);
			float decrease_0 = initial_0 * local_delta / xmax; // or xmid?
			point(0) -= decrease_0;

			float initial_1 = point(1);
			if(initial_1 > ymid) {
				float decrease_1 = (local_delta/(ymax-ymid))* (initial_1 - ymid);
				point(1) -= decrease_1;
			} else {
				float increase_1 = (local_delta/(ymid-ymin))* (ymid - initial_1);
				point(1) += increase_1;
			}


		} else if (point(2)>z3){ // section 3

			float local_delta = (delta_max / (z3-z4)) * point(2) + ((z4*delta_max)/(z4-z3));

			float initial_0 = point(0);
			float decrease_0 = initial_0 * local_delta / xmax; // or xmid?
			point(0) -= decrease_0;

			float initial_1 = point(1);
			if(initial_1 > ymid) {
				float decrease_1 = (local_delta/(ymax-ymid))* (initial_1 - ymid);
				point(1) -= decrease_1;
			} else {
				float increase_1 = (local_delta/(ymid-ymin))* (ymid - initial_1);
				point(1) += increase_1;
			}

		}
	} // ELSE DO NOTHING
}

void GridTransformation::setWrinklesParameter(Parameters& param) {
	defect_location = param.wrinklePos;
	wrinkleOri = param.wrinkleOri;

	defect_size = param.wrinkleSize;
	ref_lenght = defect_size+5.*delta_max;

	damping = {param.wrinkleDamp(0)/ref_lenght, param.wrinkleDamp(1)/(ymax-ymin), param.wrinkleDamp(2)/(zmax-zmin)}; // the more damping is big the small is the wrinkle in that direction
}

void GridTransformation::Y_alligned_wrinkles(Vector3f& point){

	Vector3f init = point;

	/* Define defect function */
	// if(init[2] < defect_location[2]-2.0*defect_size || \
	//    init[2] > defect_location[2]+2.0*defect_size ) {   //Place defect only in a section of the part
	// 	return;
	// }

	point[0] += defect_size \
		* 1./pow(cosh(damping[0] * M_PI * (init[0] - defect_location[0])),2.)\
		* 1./pow(cosh(damping[1] * M_PI * (init[1] - defect_location[1])),2.)\
		* 1./pow(cosh(damping[2] * M_PI * (init[2] - defect_location[2])),2.);


	return;
}

void GridTransformation::wrinkles(Vector3f& point){

	// Wrinkle in plan XY
	Vector3f init = point;

	float xref = defect_location[0];
	float yref = defect_location[1];
	float zref = defect_location[2] + (ymid - init[1])*tan(wrinkleOri*M_PI/180.0);

	point[0] += defect_size \
		* 1./pow(cosh(damping[0] * M_PI * (init[0] - xref)),2.)\
		* 1./pow(cosh(damping[1] * M_PI * (init[1] - yref)),2.)\
		* 1./pow(cosh(damping[2] * M_PI * (init[2] - zref)),2.);

	point[1] += 0.5 * defect_size \
		* 1./pow(cosh(damping[0] * M_PI * (init[0] - xref)),2.)\
		* 1./pow(cosh(damping[1] * M_PI * (init[1] - yref)),2.)\
		* 1./pow(cosh(damping[2] * M_PI * (init[2] - zref)),2.);


	return;
}


/////////////////////////////////////////////////////////////////////////////////////


Vector3f rotate_hors_axes(const Vector3f& u, float theta, const Vector3f& v) { // ref vector, angle, vector to be rot
	Vector3f urot;
	if (std::abs(theta)>0.0) {
		Matrix3f Q,R,id,Q2;
		Q(0,0)=0.0;
		Q(1,1)=0.0;
		Q(2,2)=0.0;
		Q(1,0)=-u(2);
		Q(2,0)=u(1);
		Q(2,1)=-u(0);
		Q(0,1)=-Q(1,0);
		Q(0,2)=-Q(2,0);
		Q(1,2)=-Q(2,1);
		id.setIdentity();
		Q2 = Q*Q;
		theta=theta*M_PI/180;
		R = id + std::sin(theta) * Q + (1-std::cos(theta)) * Q2;
		urot = R*v;
		urot.normalize();
		for (int i = 0; i < 3; i++) {
			if (std::abs(urot(i))<0.000001) {
				urot(i) = 0.0;
			}
		}
	} else {
		urot=v;
	}
	return urot;
}

void find_U(Mesh& m, Vector3f& u, Vector3f bary) {
	Vector3f Z (0.0, 0.0, 1.0);
	Vector3f Z_ (0.0, 0.0, -1.0);

	// float mid = (P_[1](0)+P_[0](0))/2.0;
	float xmin = std::min(m.P[0](0), m.P[1](0));
	float xmax = std::max(m.P[0](0), m.P[1](0));

	float ymin = std::min(m.P[0](1), m.P[1](1));
	float ymax = std::max(m.P[0](1), m.P[1](1));

	if (bary(0) < xmin) {
		if (bary(1) < ymin) {
			u(0) = 1.0;
		} else if (bary(1) > ymax) {
			u(0) = -1.0;
		}
	} else {
		if (ymin < bary(1) && bary(1) < ymax) {
			u(1) = 1.0;
		} else if (bary(1) < ymin) {
			Vector3f tmp (0.0, 0.0, 0.0);
			// tmp(0) = std::abs(bary(0) - xmin);
			// tmp(1) = std::abs(bary(1) - ymin);
			tmp(0) = bary(0) - xmin;
			tmp(1) = bary(1) - ymin;
			tmp.normalize();
			u=tmp.cross(Z_);
			// std::cout << u(0) << ", " << u(1) << ", " << u(2) << std::endl;
		} else if (bary(1) > ymax) {
			Vector3f tmp (0.0, 0.0, 0.0);
			tmp(0) = bary(0) - xmax;
			tmp(1) = bary(1) - ymax;
			tmp.normalize();
			u=tmp.cross(Z_);
		}
	}
	return;
}

void localCoorSyst(Mesh& m) {

	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			Vector3f u, v, urot, w;
			u=Vector3f::Zero();
			v=Vector3f::Zero();
			w=Vector3f::Zero();
			urot=Vector3f::Zero();

			// ~~~~~~COMPTAGES~~~~~~
			Vector3f bary = Vector3f::Zero();

			for(int i=0; i<elem.size; i++)
				bary += m.vertices.col(elem.Nodes(i,l));
			bary *= 1.0 / (elem.size+0.0);

			find_U(m, u, bary);

			Vector3f Z (0.0, 0.0, 1.0);

			// u.normalize();
			v = u.cross(Z);
			v.normalize();
			w=u.cross(v);
			w.normalize();

			for(int i=0; i<3;++i) {
				elem.U(i,l) = u(i);
				elem.V(i,l) = v(i);
				elem.W(i,l) = w(i);
			}
		}
	}
}

void GeometricTransformation(Mesh& m, Parameters& param) {

	if(param.add_ramp==false && param.add_wrinkles==false){
		std::cout << "GeoTransformation is used for nothing." << std::endl;
		return;
	}

	GridTransformation GT;
	GT.initialise(m, param);
	GT.setWrinklesParameter(param);

	for(int node=0; node < m.Nb_vertices(); ++node) {

		Vector3f point;
		point(0) = m.vertices(0, node);
		point(1) = m.vertices(1, node);
		point(2) = m.vertices(2, node);

		if(param.add_ramp)
			GT.ramp(point);
		if(param.add_wrinkles)
			GT.wrinkles(point);

		m.vertices(0, node) = point(0);
		m.vertices(1, node) = point(1);
		m.vertices(2, node) = point(2);
	}


	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			Vector3f u, v, w, urot, vrot, wrot;
			for(int i=0; i<3;++i) {
				u(i)=elem.U(i,l);
				v(i)=elem.V(i,l);
				w(i)=elem.W(i,l);
			}
			urot=Vector3f::Zero();
			vrot=Vector3f::Zero();
			wrot=Vector3f::Zero();
			// ~~~~~~COMPTAGES~~~~~~
			Vector3f bary=Vector3f::Zero();
			for(int i=0; i<elem.size; i++)
				bary += m.vertices.col(elem.Nodes(i,l));
			bary *= 1.0 / (elem.size+0.0);

			float h = 0.1;
			///
			Vector3f baryPu = bary + h * u;
			Vector3f baryMu = bary - h * u;
			if(param.add_ramp){
				GT.ramp(baryPu);
				GT.ramp(baryMu);
			}
			if(param.add_wrinkles){
				GT.wrinkles(baryPu);
				GT.wrinkles(baryMu);
			}
			Vector3f local_dir = baryPu-baryMu;
			auto angle1 = -atan2(local_dir.dot(u)/h,local_dir.dot(w)/h)*180.0/M_PI+90.0; // atan(dir0/dir2) --> adjacent dir is 2

			///
			Vector3f baryPv = bary + h * v;
			Vector3f baryMv = bary - h * v;
			if(param.add_ramp){
				GT.ramp(baryPv);
				GT.ramp(baryMv);
			}
			if(param.add_wrinkles){
				GT.wrinkles(baryPv);
				GT.wrinkles(baryMv);
			}
			local_dir = baryPv-baryMv;
			auto angle2 = -atan2(local_dir.dot(v)/h,local_dir.dot(u)/h)*180.0/M_PI+90; // atan(dir1/dir0) --> adjacent dir is 0

			///
			Vector3f baryPw = bary + h * w;
			Vector3f baryMw = bary - h * w;
			if(param.add_ramp){
				GT.ramp(baryPw);
				GT.ramp(baryMw);
			}
			if(param.add_wrinkles){
				GT.wrinkles(baryPw);
				GT.wrinkles(baryMw);
			}
			local_dir = baryPw-baryMw;
			auto angle0 = -atan2(local_dir.dot(w)/h,local_dir.dot(v)/h)*180.0/M_PI+90.0; // atan(dir2/dir1) --> adjacent dir is 1

			// if (angle2 > 0.1)
			// 	std::cout << "angle0 : " << angle0 << ", angle1 : " << angle1 << ", angle2 : " << angle2 << std::endl;

			urot = rotate_hors_axes(v, angle1, u);
			urot = rotate_hors_axes(w, angle2, urot);

			vrot = rotate_hors_axes(u, angle0, v);
			vrot = rotate_hors_axes(w, angle2, vrot);

			wrot = rotate_hors_axes(u, angle0, w);
			wrot = rotate_hors_axes(v, angle1, wrot);

			for(int i=0; i<3;++i) {
				elem.U(i,l) = urot(i);
				elem.V(i,l) = vrot(i);
				elem.W(i,l) = wrot(i);
			}
		}
	}

	return;
}

void StackingSequence(Mesh& m) {

	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			Vector3f u, v, w, urot, wrot;
			for(int i=0; i<3;++i) {
				u(i)=elem.U(i,l);
				v(i)=elem.V(i,l);
				w(i)=elem.W(i,l);
			}
			urot=Vector3f::Zero();
			wrot=Vector3f::Zero();

			Vector3f Z (0.0, 0.0, 1.0);

			urot = rotate_hors_axes(v, m.stacking_sequence[elem.Markers(0,l)-1], u);

			wrot=urot.cross(v);
			wrot.normalize();

			for(int i=0; i<3;++i) {
				elem.U(i,l) = urot(i);
				elem.W(i,l) = wrot(i);
			}
		}
	}
}