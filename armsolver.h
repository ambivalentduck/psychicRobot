#ifndef ARMSOLVER_H
#define ARMSOLVER_H

#include <gsl/gsl_errno.h>
//#define NEWGSL
#ifdef NEWGSL
	#include <gsl/gsl_odeiv2.h>
#else
	#include <gsl/gsl_odeiv.h>
#endif
#include <deque>
#include "arm.h"
#include "displaywidget.h"
#include <QSemaphore>
#include <QMutex>
#include <QThread>


class ArmSolver : public QThread
{
public:
	ArmSolver(twoLinkArm::ArmParams P, bool solveIntent=true, bool constImpedance=true);
	~ArmSolver();
	int func(double t, const double y[], double f[]);
	void setParams(twoLinkArm::ArmParams P);
	void push(double t, point p, point v, point a, point force, mat2 kp=mat2(0,0,0,0), mat2 kd=mat2(0,0,0,0));
	bool pull(point &p, int timeout=-1);
	void solve();
	void run();
	static int statfunc(double t, const double y[], double f[], void *params);
	static int statjac(double t, const double y[], double *dfdy, double dfdt[], void *params);
private:
	#ifdef NEWGSL
	gsl_odeiv2_driver * driver;
	gsl_odeiv2_system sys;
	#else
	const gsl_odeiv_step_type * odetype;
	gsl_odeiv_step * odestep;
	gsl_odeiv_control * odecontrol;
	gsl_odeiv_evolve * odeevolve;
	gsl_odeiv_system sys;
	#endif
		
	void * voidpointer;
	twoLinkArm * arm;
	bool solveDes, constImp, seeded;
	std::deque<point> qm,qs,qmdot,qsdot,qmddot,torquem;
	std::deque<mat2> Kpm,Kdm;
	std::deque<double> times, stimes;
	mat2 Kd, Kp;
	QSemaphore solvesemaphore, grabsemaphore;
	QMutex destructomutex, paramsMutex;
	point qst, qstdot;
};



#endif

