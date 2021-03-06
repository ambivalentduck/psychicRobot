#include "arm.h"
#include <cmath>
#include <iostream>

twoLinkArm::twoLinkArm(ArmParams P)
{
	params=P;
}

twoLinkArm::~twoLinkArm()
{
	return;
}

void twoLinkArm::setParams(ArmParams P)
{
	params=P;
}

void twoLinkArm::setShoulder(point x0)
{
	params.x0=x0;
}

void twoLinkArm::moveShoulder(point x0)
{
	params.x0=params.x0+x0;
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
	x=x-params.x0;
	double S=(std::pow(params.l1+params.l2,2)-(std::pow(x.X(),2)+std::pow(x.Y(),2)))/(std::pow(x.X(),2)+std::pow(x.Y(),2)-std::pow(params.l1-params.l2,2));
	if(S<0) {q=point(); return false;}
	q[1]=-2.0*std::atan(std::sqrt(S)); //Negative sign indicates elbow right.
	q[0]=std::atan2(x.Y(),x.X())-std::atan2(params.l2*std::sin(q[1]),params.l1+params.l2*std::cos(q[1]));
	return true;
}

void twoLinkArm::ikinAbs(point x, point &q)
{
	x=x-params.x0;
	double S=(std::pow(params.l1+params.l2,2)-(std::pow(x.X(),2)+std::pow(x.Y(),2)))/(std::pow(x.X(),2)+std::pow(x.Y(),2)-std::pow(params.l1-params.l2,2));
	q[1]=-2.0*std::atan(std::sqrt(std::abs(S))); //Negative sign indicates elbow right.
	q[0]=std::atan2(x.Y(),x.X())-std::atan2(params.l2*std::sin(q[1]),params.l1+params.l2*std::cos(q[1]));
}

void twoLinkArm::computeInertiaCoriolis(point q, point qdot, mat2 &D, point &C)
{
	double c2=std::cos(q.Y());
	double d11=params.m1*std::pow(params.lc1,2)+params.m2*(std::pow(params.l1,2)+std::pow(params.lc2,2)+2.0*params.l1*params.lc2*c2)+params.I1+params.I2;
	double d12=params.m2*(std::pow(params.lc2,2)+params.l1*params.lc2*c2)+params.I2;
	double d21=d12;
	double d22=params.m2*std::pow(params.lc2,2)+params.I2;
	double h=params.m2*params.l1*params.lc2*std::sin(q.Y());
	
	D=mat2(d11,d12,d21,d22);
	C=point(h*qdot.Y()*(2*qdot.X()+qdot.Y()), h*std::pow(qdot.X(),2));
}

point twoLinkArm::computeDynamics(point qDes, point qDesDot, point qDesDDot, point q, point qDot, point torque, point torqueFB)
{
	mat2 Dreal, Dexp;
	point Creal, Cexp;
	computeInertiaCoriolis(q,qDot,Dreal,Creal);
	computeInertiaCoriolis(qDes,qDesDot,Dexp,Cexp);
	point torqueFF=Dexp*qDesDDot+Cexp;
	return point(Dreal/(torqueFF-torqueFB-torque-Creal));
}

point twoLinkArm::computeInvDynamics(point q, point qDot, point qDDot, point qDes, point qDesDot, point torque, point torqueFB)
{
	mat2 Dreal, Dexp;
	point Creal, Cexp;
	computeInertiaCoriolis(q,qDot,Dreal,Creal);
	computeInertiaCoriolis(qDes,qDesDot,Dexp,Cexp);
	point torqueFFreal=Dreal*qDDot+Creal;
	return point(Dexp/(torqueFFreal+torqueFB+torque-Cexp));
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
	
	return 1;
}

twoLinkArm::ArmParams twoLinkArm::calcParams(double weight, double height, point x0)
{
	double h=0.0254*height;
	return calcParams(weight,.186*h,.1964*h,x0);
}

twoLinkArm::ArmParams twoLinkArm::calcParams(double weight, double l1, double l2, point x0)
{
	ArmParams P;
	double m=0.453592*weight;
	P.x0=x0;
	P.l1=l1;
	P.l2=l2;
	P.lc1=.436*l1;
	P.lc2=.682*l2;
	P.m1=.028*m;
	P.m2=.022*m;
	P.I1=P.m1*pow(.322*l1,2);
	P.I2=P.m2*pow(.468*l2,2);
	return P;
}



