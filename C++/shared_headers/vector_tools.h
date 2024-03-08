#pragma once

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/LU>
#include <Eigen/Geometry>

#define PI 3.14159265

using namespace Eigen;

Vector3d rotate_hors_axes(const Vector3d& axe_rot, double theta, const Vector3d& to_be_rotated) { // ref vector, angle, vector to be rot
	Vector3d urot;
	if (std::abs(theta)>0.0) {
		Matrix3d Q,R,id,Q2;
		Q(0,0)=0.0;
		Q(1,1)=0.0;
		Q(2,2)=0.0;
		Q(1,0)=-axe_rot(2);
		Q(2,0)= axe_rot(1);
		Q(2,1)=-axe_rot(0);
		Q(0,1)=-Q(1,0);
		Q(0,2)=-Q(2,0);
		Q(1,2)=-Q(2,1);
		id.setIdentity();
		Q2 = Q*Q;
		theta=theta*PI/180;
		R = id + std::sin(theta) * Q + (1-std::cos(theta)) * Q2;
		urot = R*to_be_rotated;
		urot.normalize();
		for (int i = 0; i < 3; i++) {
			if (std::abs(urot(i))<0.000001) {
				urot(i) = 0.0;
			}
		}
	} else {
		urot=to_be_rotated;
	}
	return urot;
}
