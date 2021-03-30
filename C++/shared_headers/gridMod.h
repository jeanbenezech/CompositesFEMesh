#pragma once

#include "mesh.h"
#include "element.h"
#include "parameters.h"
#include "vector_tools.h"

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/LU>
#include <Eigen/Geometry>

#define PI 3.14159265

using namespace Eigen;

class GridTransformation{
	public:

	// Global gemoetry
	double z1, z2, z3, z4;
	double delta_max;
	double zmin, zmax, ymin, ymax, xmin, xmax;
	double ymid;

	// Wrinkles parameters
	Vector3d defect_location, damping, center_rot;
	double wrinkleOri;
	double ref_lenght;
	double defect_size;

	// General parameters
	int verbosity = 0;

	void initialise(Mesh& m, Parameters& param);
	void ramp(Vector3d& point);
	void setWrinklesParameter(Parameters& param);
	void OneDir_wrinkles(Vector3d& point, int dim);
	void wrinkles(Vector3d& point);

	private:
	// double epsilon = 0.1;
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

	double deltaz = (zmax-zmin) / 16.0;
	z1 = zmin + 2.0*deltaz;
	z2 = zmin + 7.0*deltaz;
	z3 = zmin + 9.0*deltaz;
	z4 = zmin + 14.0*deltaz;

	delta_max = param.rampSize;

	ymid = (ymax-ymin)/2.0;

	// std::cout << "Ymid : " << ymid << std::endl;
}

void GridTransformation::ramp(Vector3d& point){
	double ymid = (ymax+ymin)/2.0;

	if (point(2)>z1 && point(2)<z4){ // section 1, 2, 3

		if (point(2)>=z2 && point(2)<=z3){ // section 2

			double initial_0 = point(0);
			double decrease_0 = initial_0 * delta_max / xmax; // or xmid?
			point(0) -= decrease_0;

			double initial_1 = point(1);
			if(initial_1 > ymid) {
				double decrease_1 = (delta_max/(ymax-ymid))* (initial_1 - ymid);
				point(1) -= decrease_1;
			} else {
				double increase_1 = (delta_max/(ymid-ymin))* (ymid - initial_1);
				point(1) += increase_1;
			}

		} else if (point(2)<z2){ // section 1

			double local_delta = (delta_max / (z2-z1)) * point(2) + ((z1*delta_max)/(z1-z2));

			double initial_0 = point(0);
			double decrease_0 = initial_0 * local_delta / xmax; // or xmid?
			point(0) -= decrease_0;

			double initial_1 = point(1);
			if(initial_1 > ymid) {
				double decrease_1 = (local_delta/(ymax-ymid))* (initial_1 - ymid);
				point(1) -= decrease_1;
			} else {
				double increase_1 = (local_delta/(ymid-ymin))* (ymid - initial_1);
				point(1) += increase_1;
			}


		} else if (point(2)>z3){ // section 3

			double local_delta = (delta_max / (z3-z4)) * point(2) + ((z4*delta_max)/(z4-z3));

			double initial_0 = point(0);
			double decrease_0 = initial_0 * local_delta / xmax; // or xmid?
			point(0) -= decrease_0;

			double initial_1 = point(1);
			if(initial_1 > ymid) {
				double decrease_1 = (local_delta/(ymax-ymid))* (initial_1 - ymid);
				point(1) -= decrease_1;
			} else {
				double increase_1 = (local_delta/(ymid-ymin))* (ymid - initial_1);
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

void GridTransformation::OneDir_wrinkles(Vector3d& point, int dim){

	Vector3d init = point;

	double xref = defect_location[0];
	double yref = defect_location[1];
	double zref = defect_location[2] + (ymid - init[1])*tan(wrinkleOri*PI/180.0);

	point[dim] += defect_size \
		* 1./pow(cosh(damping[0] * PI * (init[0] - xref)),2.)\
		* 1./pow(cosh(damping[1] * PI * (init[1] - yref)),2.)\
		* 1./pow(cosh(damping[2] * PI * (init[2] - zref)),2.);

	return;
}

void GridTransformation::wrinkles(Vector3d& point){

	// Wrinkle in plan XY
	Vector3d init = point;

	double xref = defect_location[0];
	double yref = defect_location[1];
	double zref = defect_location[2] + (ymid - init[1])*tan(wrinkleOri*PI/180.0);

	point[0] += defect_size \
		* 1./pow(cosh(damping[0] * PI * (init[0] - xref)),2.)\
		* 1./pow(cosh(damping[1] * PI * (init[1] - yref)),2.)\
		* 1./pow(cosh(damping[2] * PI * (init[2] - zref)),2.);

	point[1] += 0.5 * defect_size \
		* 1./pow(cosh(damping[0] * PI * (init[0] - xref)),2.)\
		* 1./pow(cosh(damping[1] * PI * (init[1] - yref)),2.)\
		* 1./pow(cosh(damping[2] * PI * (init[2] - zref)),2.);


	return;
}

/////////////////////////////////////////////////////////////////////////////////////

void find_U(Mesh& m, Vector3d& u, Vector3d bary) {
	Vector3d Z (0.0, 0.0, 1.0);
	Vector3d Z_ (0.0, 0.0, -1.0);

	// double mid = (P_[1](0)+P_[0](0))/2.0;
	double xmin = std::min(m.P[0](0), m.P[1](0));
	double xmax = std::max(m.P[0](0), m.P[1](0));

	double ymin = std::min(m.P[0](1), m.P[1](1));
	double ymax = std::max(m.P[0](1), m.P[1](1));

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
			Vector3d tmp (0.0, 0.0, 0.0);
			// tmp(0) = std::abs(bary(0) - xmin);
			// tmp(1) = std::abs(bary(1) - ymin);
			tmp(0) = bary(0) - xmin;
			tmp(1) = bary(1) - ymin;
			tmp.normalize();
			u=tmp.cross(Z_);
			// std::cout << u(0) << ", " << u(1) << ", " << u(2) << std::endl;
		} else if (bary(1) > ymax) {
			Vector3d tmp (0.0, 0.0, 0.0);
			tmp(0) = bary(0) - xmax;
			tmp(1) = bary(1) - ymax;
			tmp.normalize();
			u=tmp.cross(Z_);
		}
	}
	return;
}

void localCoorSyst(Mesh& m, Parameters& param) {

	for(auto& elem : m.Elements){ // Loop over element types
		for(int l=0; l < elem.nb; ++l) { // Loop over elements for each type

			Vector3d u, v, urot, w;
			u=Vector3d::Zero();
			v=Vector3d::Zero();
			w=Vector3d::Zero();
			urot=Vector3d::Zero();

			if(elem.Markers(0,l) == param.resin_id){ // Resin layer got global coord sys
				elem.U(0,l) = 1.0;
				elem.U(1,l) = 0.0;
				elem.U(2,l) = 0.0;

				elem.V(0,l) = 0.0;
				elem.V(1,l) = 1.0;
				elem.V(2,l) = 0.0;

				elem.W(0,l) = 0.0;
				elem.W(1,l) = 0.0;
				elem.W(2,l) = 1.0;

			} else {
				Vector3d bary = Vector3d::Zero();
				for(int i=0; i<elem.size; i++)
					bary += m.vertices.col(elem.Nodes(i,l));
				bary *= 1.0 / (elem.size+0.0);

				find_U(m, u, bary);

				Vector3d Z (0.0, 0.0, 1.0);

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
}

void globalCoorSyst(Mesh& m) {

	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			elem.U(0,l) = 1.0;
			elem.U(1,l) = 0.0;
			elem.U(2,l) = 0.0;

			elem.V(0,l) = 0.0;
			elem.V(1,l) = 1.0;
			elem.V(2,l) = 0.0;

			elem.W(0,l) = 0.0;
			elem.W(1,l) = 0.0;
			elem.W(2,l) = 1.0;

		}
	}
}

void GeometricTransformation(Mesh& m, Parameters& param) {

	if(param.add_ramp==false && param.add_wrinkles==false){
		std::cout << "GeoTransformation is used for nothing." << std::endl;
		return;
	}

	int dim = 0;
	if(param.Shape==1)
		dim = 1;

	GridTransformation GT;
	GT.initialise(m, param);
	GT.setWrinklesParameter(param);

	for(int node=0; node < m.Nb_vertices(); ++node) {

		Vector3d point;
		point(0) = m.vertices(0, node);
		point(1) = m.vertices(1, node);
		point(2) = m.vertices(2, node);

		if(param.add_ramp)
			GT.ramp(point);
		if(param.add_wrinkles)
			GT.OneDir_wrinkles(point, dim);


		m.vertices(0, node) = point(0);
		m.vertices(1, node) = point(1);
		m.vertices(2, node) = point(2);
	}


	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			Vector3d u, v, w, urot, vrot, wrot;
			for(int i=0; i<3;++i) {
				u(i)=elem.U(i,l);
				v(i)=elem.V(i,l);
				w(i)=elem.W(i,l);
			}
			urot=Vector3d::Zero();
			vrot=Vector3d::Zero();
			wrot=Vector3d::Zero();
			// ~~~~~~COMPTAGES~~~~~~
			Vector3d bary=Vector3d::Zero();
			for(int i=0; i<elem.size; i++)
				bary += m.vertices.col(elem.Nodes(i,l));
			bary *= 1.0 / (elem.size+0.0);

			double h = 0.1;

			///
			Vector3d baryPu = bary + h * u;
			Vector3d baryMu = bary - h * u;
			if(param.add_ramp){
				GT.ramp(baryPu);
				GT.ramp(baryMu);
			}
			if(param.add_wrinkles){
				GT.OneDir_wrinkles(baryPu, dim);
				GT.OneDir_wrinkles(baryMu, dim);
			}

			///
			Vector3d baryPv = bary + h * v;
			Vector3d baryMv = bary - h * v;
			if(param.add_ramp){
				GT.ramp(baryPv);
				GT.ramp(baryMv);
			}
			if(param.add_wrinkles){
				GT.OneDir_wrinkles(baryPv, dim);
				GT.OneDir_wrinkles(baryMv, dim);
			}

			///
			Vector3d baryPw = bary + h * w;
			Vector3d baryMw = bary - h * w;
			if(param.add_ramp){
				GT.ramp(baryPw);
				GT.ramp(baryMw);
			}
			if(param.add_wrinkles){
				GT.OneDir_wrinkles(baryPw, dim);
				GT.OneDir_wrinkles(baryMw, dim);
			}


			double angle0=0.0, angle1=0.0, angle2=0.0;


			// No shift in the V direction which is the normal to the ply
			// if(dim==1){
			Vector3d local_dir_use_U = baryPu-baryMu;
			angle1 = -atan2(local_dir_use_U.dot(u)/h,local_dir_use_U.dot(w)/h)*180.0/PI+90.0; // atan(dir0/dir2) --> adjacent dir is 2

			// Vector3d local_dir_use_V = baryPv-baryMv;
			angle2 = -atan2(local_dir_use_U.dot(v)/h,local_dir_use_U.dot(u)/h)*180.0/PI;

			Vector3d local_dir_use_W = baryPw-baryMw;
			angle0 = -atan2(local_dir_use_W.dot(w)/h,local_dir_use_W.dot(v)/h)*180.0/PI+90.0; // atan(dir2/dir1) --> adjacent dir is 1

			// ///
			// Vector3d local_dir_use_U = baryPu-baryMu;
			// angle1 = -atan2(local_dir_use_U.dot(u)/h,local_dir_use_U.dot(w)/h)*180.0/PI+90.0; // atan(dir0/dir2) --> adjacent dir is 2

			// ///
			// Vector3d local_dir_use_V = baryPv-baryMv;
			// angle2 = -atan2(local_dir_use_V.dot(v)/h,local_dir_use_V.dot(u)/h)*180.0/PI+90.0; // atan(dir1/dir0) --> adjacent dir is 0

			// ///
			// Vector3d local_dir_use_W = baryPw-baryMw;
			// angle0 = -atan2(local_dir_use_W.dot(w)/h,local_dir_use_W.dot(v)/h)*180.0/PI+90.0; // atan(dir2/dir1) --> adjacent dir is 1

			// if (angle1 > 0.2)
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

void StackingSequence(Mesh& m, Parameters& param) {

	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			if(elem.Markers(0,l) != param.resin_id && elem.Markers(0,l) != param.cz_id){ // Resin layer or CZ do not need StaqSeq
				Vector3d u, v, w, urot, wrot;
				for(int i=0; i<3;++i) {
					u(i)=elem.U(i,l);
					v(i)=elem.V(i,l);
					w(i)=elem.W(i,l);
				}
				urot=Vector3d::Zero();
				wrot=Vector3d::Zero();

				Vector3d Z (0.0, 0.0, 1.0);

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
}