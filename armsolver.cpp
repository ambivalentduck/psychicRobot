#include "armsolver.h"

ArmSolver::ArmSolver(twoLinkArm::ArmParams P, SOLVETYPE SolType)
{
	solType=SolType;
	impSeeded=false;
	arm=new twoLinkArm(P);
	voidpointer=(void*) this;
	sys = {statfunc, statjac, 4, voidpointer};
	#ifdef NEWGSL
	driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,1e-6, 1e-6, 0.0);
	#else
	odetype=gsl_odeiv_step_rk8pd;
	odestep=gsl_odeiv_step_alloc (odetype, 4);
	odecontrol=gsl_odeiv_control_y_new (1e-6, 0.0);
	odeevolve=gsl_odeiv_evolve_alloc(4);
	#endif
	seeded=false;
}

ArmSolver::~ArmSolver()
{
	#ifdef NEWGSL
	gsl_odeiv2_driver_free(driver);
	#else
	gsl_odeiv_evolve_free(odeevolve);
	gsl_odeiv_control_free(odecontrol);
	gsl_odeiv_step_free(odestep);
	#endif
	delete arm;
}

int ArmSolver::func(double t, const double y[], double f[])
{
	double s=(t-times[0])/(times[1]-times[0]);
	double oms=1l-s;
	point qmi=qm[0]*oms+qm[1]*s;
	point qmdoti=qmdot[0]*oms+qmdot[1]*s;
	point qmddoti=qmddot[0]*oms+qmddot[1]*s;
	point torquemi=torquem[0]*oms+torquem[1]*s;
	
	mat2 kpi;
	mat2 kdi;
	if(constImp) {kpi=Kp; kdi=Kd;}
	else {kpi=Kpm[0]*oms+Kpm[1]*s; kdi=Kdm[0]*oms+Kdm[1]*s;}
	
	f[0]=y[2];
	f[1]=y[3];
	
	point qsDDot;
	
	paramsMutex.lock();
	if(solveDes) qsDDot=arm->computeInvDynamics(qmi, qmdoti, qmddoti, point(y[0],y[1]), point(y[2],y[3]), torquemi, kpi, kdi);
	else qsDDot=arm->computeDynamics(qmi, qmdoti, qmddoti, point(y[0],y[1]), point(y[2],y[3]), torquemi, kpi, kdi);
	paramsMutex.unlock();
	
	f[2]=qsDDot[0];
	f[3]=qsDDot[1];
	
	return GSL_SUCCESS;
}

void ArmSolver::setParams(twoLinkArm::ArmParams P)
{	
	paramsMutex.lock();
	arm->setParams(P); 
	paramsMutex.unlock();
}

void ArmSolver::push(double t, point p, point v, point a, point force, mat2 kp, mat2 kd)
{
	destructomutex.lock(); //We block here just in case we push at the moment of a crash.
	point q;
	if(arm->ikin(p,q)) questionable.push_back(false);
	else {arm->ikinAbs(p,q); questionable.push_back(true);}
	qm.push_back(q);
	
	mat2 fJ=arm->jacobian(q);
	point qdot=fJ/v;
	qmdot.push_back(qdot);
	
	point qddot=arm->getQDDot(q,qdot,a);
	qmddot.push_back(qddot);
	
	point torque=(fJ.transpose())*force;
	torquem.push_back(torque);
	
	times.push_back(t);
	
	if(!constImp)
	{
		Kpm.push_back(kp);
		Kdm.push_back(kd);
	}
	
	if(!seeded)
	{
		qst=q;
		qstdot=qdot;
		if(constImp){Kp=kp; Kd=kd;}
		seeded=true;
	}
	else solvesemaphore.release();	
	destructomutex.unlock();
}

void ArmSolver::solve()
{
	if(!isRunning())
		start();
}

void ArmSolver::run()
{
	while(true)
	{
		if(!solvesemaphore.tryAcquire(1,17)) break; //Try for at least 1 frame to start solving
		double y[4];
		y[0]=qst[0];
		y[1]=qst[1];
		y[2]=qstdot[0];
		y[3]=qstdot[1];
		double t=times[0];
		double tk=times[1];
		#ifdef NEWGSL
		int status=gsl_odeiv2_driver_apply(driver, &t, tk, y);
		#else
		double h=1e-6;
		int status;
		while(t<tk)
		{
			status = gsl_odeiv_evolve_apply(odeevolve,odecontrol,odestep,&sys,&t,tk,&h,y);
			if(status!=GSL_SUCCESS) break; //Don't lock up on a crash
		}
		#endif
		if(status!=GSL_SUCCESS)
		{
			std::cout<<"Oops. Cuh-rash."<<std::endl;
			destructomutex.lock(); //We lock a mutex so that we can safely remove extant solutions.
			qm.clear();
			qmdot.clear();
			qmddot.clear();
			torquem.clear();
			times.clear();
			if(!constImp)
			{
				Kpm.clear();
				Kdm.clear();
			}
			seeded=false;
			//Since we just emptied the queues, also empty the semaphore
			solvesemaphore.acquire(solvesemaphore.available());
			destructomutex.unlock();
			break; //We just emptied the semaphore
		}
		
		qst=point(y[0],y[1]);
		qstdot=point(y[2],y[3]);
		
		qs.push_back(qst);
		qsdot.push_back(qstdot);
		stimes.push_back(t);
		questionableGrab.push_back(questionable.front());
		grabsemaphore.release();
		
		qm.pop_front();
		qmdot.pop_front();
		qmddot.pop_front();
		torquem.pop_front();
		times.pop_front();
		questionable.pop_front();
		if(!constImp)
		{
			Kpm.pop_front();
			Kdm.pop_front();
		}
	}
}

bool ArmSolver::pull(point &p, bool &dodgy, int timeout)
{
	if(!grabsemaphore.tryAcquire(1,timeout)) return false;
	p=arm->fkin(qs.front()); //Update with braindead q+qdot*delta(t)?
	dodgy=questionableGrab.front();
	questionableGrab.pop_front();
	qs.pop_front();
	qsdot.pop_front();
	stimes.pop_front();
	return true;
}

int ArmSolver::statjac(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	return GSL_SUCCESS;
}

int ArmSolver::statfunc(double t, const double y[], double f[], void *params)
{
	ArmSolver * o=(ArmSolver*) params;
	return o->func(t,y,f);
}
		
