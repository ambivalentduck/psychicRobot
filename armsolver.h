#ifndef ARMSOLVER_H
#define ARMSOLVER_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <deque>
#include "arm.h"
#include "displaywidget.h"
#include <QSemaphore>

class ArmSolver
{
public:
	ArmSolver(twoLinkArm::ArmParams P, bool solveIntent=true, bool constImpedance=true);
	~ArmSolver();
	int func(double t, const double y[], double f[]);
	void cleanpush(twoLinkArm::ArmParams P, double t, point p, point v, point a, mat2 Kp, mat2 Kd);
	void push(double t, point p, point v, point a);
	bool pull(point &p, int timeout=-1);
	void solve();
	static int statfunc(double t, const double y[], double f[], void *params);
	static int statjac(double t, const double y[], double *dfdy, double dfdt[], void *params);
private:
	gsl_odeiv2_driver * driver;
	gsl_odeiv2_system sys;
	void * voidpointer;
	twoLinkArm * arm;
	bool solveDes, constImp, seeded;
	std::deque<point> qm,qs,qmdot,qsdot,qmddot,torquem,Kpm,Kdm;
	std::deque<double> times, stimes;
	mat2 Kd, Kp;
	QSemaphore solvesemaphore, grabsemaphore;
};

class solveThread : public QThread
{
public:
	solveThread(ArmSolver * solver);
	void run();		
	ArmSolver * as;
};

#endif

