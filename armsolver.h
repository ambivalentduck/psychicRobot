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
	enum Models {CONSTIMP=0, TORQUESCALEDIMP=1, TORQUESCALEDANDREFLEX=2} model;
	
	ArmSolver(twoLinkArm::ArmParams P, Models m=TORQUESCALEDANDREFLEX, bool solveIntent=true);
	~ArmSolver();
	int func(double t, const double y[], double f[]);
	void setParams(twoLinkArm::ArmParams P);
	void setModel(Models m);
	void setImpedanceGain(double g) {impedanceGain=g;}
	void push(double t, point p, point v, point a, point force);
	bool pull(point &p, bool &dodgy, int timeout=-1);
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
	bool solveDes, seeded;
	std::deque<point> qm,qs,qmdot,qsdot,qmddot,torquem,es,esdot;
	std::deque<double> times, stimes, etimes;
	std::deque<bool> questionable, questionableGrab;
	QSemaphore solvesemaphore, grabsemaphore;
	QMutex destructomutex, paramsMutex;
	point qst, qstdot;
	double impedanceGain;
};



#endif

