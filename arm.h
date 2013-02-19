#ifndef ARM_H
#define ARM_H

#include "point.h"
#include "mat2.h"

class twoLinkArm
{
public:
	twoLinkArm(ArmParams P);
	~twoLinkArm();
	void moveShoulder(point x0);
	point fkin(point q);
	void fkin(point q, point &x, point &x1);
	bool ikin(point x, point &q);
	mat2 jacobian(point q);
	void computeInertiaCoriolis(point q, point qdot, mat2 &D, point &C);
	point computeDynamics(point qDes, point qDesDot, point qDesDDot, point q, point qDot, point torque, mat2 Kp, mat2 Kd);
	point computeInvDynamics(point q, point qDot, point qDDot, point qDes, point qDesDot, point torque, mat2 Kp, mat2 Kd);
	point getQDDot(point q, point qdot, point xddot);
	void crapPoint(point p);
	void crapMat(mat2 m);
	bool unitTests();
	
	struct ArmParams
	{
		point x0;
		double l1;
		double l2;
		double lc1;
		double lc2;
		double m1;
		double m2;
		double I1;
		double I2;
	};
	static ArmParams defaultParams();
	
private:
	 ArmParams params;
		
};

#endif

