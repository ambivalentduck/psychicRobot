#include "arm.h"
#include <cmath>
#include <iostream>

twoLinkArm::twoLinkArm(point x0, double l1, double l2, double lc1, double lc2, double m1, double m2, double I1, double I2)
{
	params.x0=x0;
	params.l1=l1;
	params.l2=l2;
	params.lc1=lc1;
	params.lc2=lc2;
	params.m1=m1;
	params.m2=m2;
	params.I1=I1;
	params.I2=I2;
}

twoLinkArm::~twoLinkArm()
{
	return;
}

void twoLinkArm::moveShoulder(point x0)
{
	params.x0=x0;
}

point twoLinkArm::fkin(point q)
{
	return point(params.x0.X()+params.l1*std::cos(q.X())+params.l2*std::cos(q.X()+q.Y()),
		params.x0.Y()+params.l1*std::sin(q.X())+params.l2*std::sin(q.X()+q.Y()));
}

void twoLinkArm::fkin(point q, point &x, point &x1)
{
	x1=point(params.x0.X()+params.l1*std::cos(q.X()),params.x0.Y()+params.l1*std::sin(q.X()));
	x=x1+point(params.l2*std::cos(q.X()+q.Y()),params.l2*std::sin(q.X()+q.Y()));
}

bool twoLinkArm::ikin(point x, point &q)
{
	point xtemp=x;
	x=x-params.x0;
	double S=(std::pow(params.l1+params.l2,2)-(std::pow(x.X(),2)+std::pow(x.Y(),2)))/(std::pow(x.X(),2)+std::pow(x.Y(),2)-std::pow(params.l1-params.l2,2));
	if(S<0) {q=point(); return false;}
	q[1]=-2.0*std::atan(std::sqrt(S));
	q[0]=std::atan2(x.Y(),x.X())-std::atan2(params.l2*std::sin(q[1]),params.l1+params.l2*std::cos(q[1]));
	return true;
}

void twoLinkArm::computeInertiaCoriolis(point q, point qdot, mat2 &D, point &C)
{
	double c2=std::cos(q.Y());
	double d11=params.m1*std::pow(params.lc1,2)+params.m2*(std::pow(params.l1,2)+std::pow(params.lc2,2)+2.0*params.l1*params.lc2*c2)+params.I1+params.I2;
	double d12=params.m2*(std::pow(params.lc2,2)+params.l1*params.lc2*c2)+params.I2;
	double d21=d12;
	double d22=params.m2*std::pow(params.lc2,2)+params.I2;
	double h=-params.m2*params.l1*params.lc2*std::sin(q.Y());
	
	D=mat2(d11,d12,d21,d22);
	C=point(2.0*h*qdot.X()*qdot.Y()+h*std::pow(qdot.Y(),2), -h*std::pow(qdot.X(),2));
}

point twoLinkArm::computeDynamics(point qDes, point qDesDot, point qDesDDot, point q, point qDot, point torque, mat2 Kp, mat2 Kd)
{
	mat2 Dreal, Dexp;
	point Creal, Cexp;
	computeInertiaCoriolis(q,qDot,Dreal,Creal);
	computeInertiaCoriolis(qDes,qDesDot,Dexp,Cexp);
	point torqueFF=Dexp*qDesDDot+Cexp;
	point torqueFB=Kp*(qDes-q)+Kd*(qDesDot-qDot);
	return point(Dreal/(torqueFF+torqueFB-torque-Creal));
}

point twoLinkArm::computeInvDynamics(point q, point qDot, point qDDot, point qDes, point qDesDot, point torque, mat2 Kp, mat2 Kd)
{
	mat2 Dreal, Dexp;
	point Creal, Cexp;
	computeInertiaCoriolis(q,qDot,Dreal,Creal);
	computeInertiaCoriolis(qDes,qDesDot,Dexp,Cexp);
	point torqueFFreal=Dreal*qDDot+Creal;
	point torqueFB=Kp*(qDes-q)+Kd*(qDesDot-qDot);
	return point(Dexp/(torqueFFreal-torqueFB+torque-Cexp));
}

mat2 twoLinkArm::jacobian(point q)
{
	double s12=std::sin(q.X()+q.Y());
	double c12=std::cos(q.X()+q.Y());
	double fJ11=-params.l1*std::sin(q.X())-params.l2*s12;
	double fJ21=params.l1*std::cos(q.X())+params.l2*c12;
	double fJ12=-params.l2*s12;
	double fJ22=params.l2*c12;
	return mat2(fJ11,fJ12,fJ21,fJ22);
}

point twoLinkArm::getQDDot(point q, point qdot, point xddot)
{
	double s12=std::sin(q.X()+q.Y());
	double c12=std::cos(q.X()+q.Y());
	double fJt1dx2=-params.l1*std::cos(q.X())-params.l2*c12;
	double fJt1dy2=-params.l2*std::cos(q.X()+q.Y());
	double fJt2dx2=-params.l1*std::sin(q.X())-params.l2*s12;
	double fJt2dy2=-params.l2*s12;
	double fJt1dxdy=-params.l2*c12;
	double fJt2dxdy=-params.l2*s12;
	return jacobian(q)/(xddot-(mat2(fJt1dx2*qdot.X()+fJt1dxdy*qdot.Y(),fJt1dy2*qdot.Y()+fJt1dxdy*qdot.X(),
		fJt2dx2*qdot.X()+fJt2dxdy*qdot.Y(),fJt2dy2*qdot.Y()+fJt2dxdy*qdot.X())*qdot));
}

void twoLinkArm::crapPoint(point p)
{
	std::cout << p[0] << std::endl << p[1] << std::endl;
}

void twoLinkArm::crapMat(mat2 m)
{
	std::cout << m.m[0] << "\t" << m.m[1] << std::endl << m.m[2] << "\t" << m.m[3] << std::endl;
}

bool twoLinkArm::unitTests()
{
	point testpoint=point(-.0908,-2.2786);
	point fcheck=fkin(testpoint);
	std::cout<<"Fcheck: "<<std::endl;
	crapPoint(fcheck);
	std::cout<<std::endl;
	
	point fcheck2;
	std::cout<<"Fcheck2: "<<std::endl;
	fkin(testpoint, fcheck, fcheck2);
	crapPoint(fcheck);
	crapPoint(fcheck2);
	std::cout<<std::endl;
	
	point icheck;
	std::cout<<"Icheck: ";
	std::cout<<ikin(fcheck,icheck)<<std::endl;;
	crapPoint(icheck);
	std::cout<<std::endl;
	
	std::cout<<"Jacobian: "<<std::endl;
	mat2 testjac=jacobian(point(5,6));
	crapMat(testjac);
	std::cout<<std::endl;
	
	std::cout<<"Q double dot: "<<std::endl;
	point qdd=getQDDot(point(1,2),point(3,4),point(5,6));
	crapPoint(qdd);
	std::cout<<std::endl;
	
	std::cout<<"DC: "<<std::endl;
	mat2 D;
	point C;
	computeInertiaCoriolis(point(1,2),point(3,4),D,C);
	crapMat(D);
	crapPoint(C);
	std::cout<<std::endl;
	
	return 0;
}


