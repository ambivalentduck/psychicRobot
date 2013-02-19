#ifndef ARMSOLVER_H
#define ARMSOLVER_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <deque>
#include "arm.h"

class ArmSolver
{
public:
	ArmSolver(twoLinkArm::ArmParams P, bool solveQ=true);
	~ArmSolver();
	int func(double t, const double y[], double f[]);
	void solve(double t1, double t2, double y1, int n);
	static int statfunc(double t, const double y[], double f[], void *params);
	static int statjac(double t, const double y[], double *dfdy, double dfdt[], void *params);
	
private:
	gsl_odeiv2_driver * driver;
	gsl_odeiv2_system sys;
	std::deque<point> solved;
	std::deque<double> times;
	void * voidpointer;
	twoLinkArm * arm;
	bool doQ;
	point
};

#endif

