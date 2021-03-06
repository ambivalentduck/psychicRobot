#ifndef ARM_H
#define ARM_H

#include "point.h"
#include "mat2.h"

class twoLinkArm
{
public:
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
	
	twoLinkArm(ArmParams P);
	~twoLinkArm();
	void moveShoulder(point x0);
	void setShoulder(point x0);
	void setParams(ArmParams P);
	point fkin(point q);
	void fkin(point q, point &x, point &x1);
	bool ikin(point x, point &q);
	void ikinAbs(point x, point &q);
	mat2 jacobian(point q);
	void computeInertiaCoriolis(point q, point qdot, mat2 &D, point &C);
	point computeDynamics(point qDes, point qDesDot, point qDesDDot, point q, point qDot, point torque, point torqueFB);
	point computeInvDynamics(point q, point qDot, point qDDot, point qDes, point qDesDot, point torque, point torqueFB);
	point getQDDot(point q, point qdot, point xddot);
	void crapPoint(point p);
	void crapMat(mat2 m);
	bool unitTests();
	
	static ArmParams calcParams(double weight, double height, point x0);
	static ArmParams calcParams(double weight, double l1, double l2, point x0);
private:
	 ArmParams params;
		
};

#endif

