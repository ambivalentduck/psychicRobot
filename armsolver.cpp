#include "armsolver.h"
#include <cstdio>

#define REFLEXDELAY .05l

ArmSolver::ArmSolver(twoLinkArm::ArmParams P, Models m, bool solveIntent)
{
	model=m;
	solveDes=solveIntent;
	arm=new twoLinkArm(P); //This is where you should set the arm model, no switching later.
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
	impedanceGain=1;
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
	//Standard linear interpolation
	double s=(t-times[0])/(times[1]-times[0]);
	double oms=1l-s;
	point qmi=qm[0]*oms+qm[1]*s;
	point qmdoti=qmdot[0]*oms+qmdot[1]*s;
	point qmddoti=qmddot[0]*oms+qmddot[1]*s;
	point torquemi=torquem[0]*oms+torquem[1]*s;
	point eri, erdoti; //initialize reflex/delayed errors to zeros
	int etimesSize=etimes.size();
	
	if(model==TORQUESCALEDANDREFLEX)
	{
		/*Because sampling is potentially extremely uneven, search starting from the back until the NEXT time both exists and is larger than t-delay
		We're quietly assuming that the front end of the array is getting trimmed intelligently. Probing for t=-1 would find 0<t<1 from t=[0 1 2...] */
		bool flag=false;
		int k=0;
		while((k+2)<etimesSize) //Restrict ourselves to points we actually have. Remember that 0-indexing means k==0 implies size==2.
		{
			if((t-REFLEXDELAY)<=etimes[k+1]) //Handles the corner case t-delay==a sampling time without requiring an additional sample
			{
				flag=true;
				break;
			}
			k++;
		}
		
		if(flag) //We have two points: the upper of which is past t-delay, the lower of which is not. 
		{
			s=((t-REFLEXDELAY)-etimes[k])/(etimes[k+1]-etimes[k]);
			oms=1l-s;
			eri=es[k]*oms+es[k+1]*s;
			erdoti=esdot[k]*oms+esdot[k+1]*s;
		}
		else //Absent data, assume error was zero.
		{
			eri=point(0,0);
			erdoti=point(0,0);
		}
	}
	
	point torqueFB;
	point torqueReflex=point(0,0);
	point e, edot;
	//Direction of solve changes which variables contain real vs desired trajectories
	if(solveDes) {e=qmi-point(y[0],y[1]); edot=qmdoti-point(y[2],y[3]);}
	else {e=point(y[0],y[1])-qmi; edot=point(y[2],y[3])-qmdoti;}
	
	mat2 Dreal,kp,kp0,kp1;
	point jointTorques, Creal;
	double jt1, jt2;
	
	switch(model)
	{
	case CONSTIMP:
		torqueFB=mat2(8,2,2,5)*e+mat2(2.3, .09, .09, 2.4)*edot;
		break;
	case TORQUESCALEDANDREFLEX:
		torqueReflex=(eri+erdoti*2l)/50l;
		//Explicitly choose NOT to call break so that we fall through into the next case.
	case TORQUESCALEDIMP:
		arm->computeInertiaCoriolis(qmi,qmdoti,Dreal,Creal);
		/*Sign conventions work out so that muscle-generated torques (FF+FB) = external + inertial torques.
		In the case that you measure the physical acceleration of the arm, this is well and good.
		In the case that you know the intention a priori, you don't have access to the inertial torque OR the feedback torque until you've solved.
		We look back a timestep and make the assumption that arm torque isn't changing quickly relative to our sampling rate.
		Still, our simulation of realization from intent requires assumptions and/or hacks.*/
		if(solveDes) jointTorques=Dreal*qmddoti+Creal+torquemi;
		else jointTorques=Dreal*qstdot+Creal+torquemi;
		jt1=jointTorques[0];
		jt2=jointTorques[1];
		jointTorques=point(jt1,jt2);
		kp0=mat2(10.8,2.83,2.51,8.67);
		kp1=mat2(3.18l*jointTorques.X(),2.15l*jointTorques.Y(),2.34l*jointTorques.X(),6.18l*jointTorques.Y());
		kp=kp0+kp1;
		torqueFB=kp*(e+edot/12l+torqueReflex); //Unless we fell through, reflex torque is 0.
		break;
	}
	torqueFB=torqueFB*impedanceGain; //I think this is just multiplication by 1.
	
	f[0]=y[2];
	f[1]=y[3];
	
	point qsDDot;
	
	paramsMutex.lock();
		
	if(solveDes) qsDDot=arm->computeInvDynamics(qmi, qmdoti, qmddoti, point(y[0],y[1]), point(y[2],y[3]), torquemi, torqueFB);
	else qsDDot=arm->computeDynamics(qmi, qmdoti, qmddoti, point(y[0],y[1]), point(y[2],y[3]), torquemi, torqueFB);
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

void ArmSolver::setModel(Models m)
{	
	paramsMutex.lock();
	model=m;
	paramsMutex.unlock();
}

void ArmSolver::push(double t, point p, point v, point a, point force)
{
	destructomutex.lock(); //We block here just in case we push at the moment of a crash.
	point q;
	if(arm->ikin(p,q)) questionable.push_back(false); //Deal with what's usually trunk or shoulder motility allowing impossible reaches.
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
	
	if(!seeded) //Set the initial conditions to the measured conditions
	{
		qst=q;
		qstdot=qdot;
		seeded=true;
	}
	else solvesemaphore.release();	
	destructomutex.unlock();
}

void ArmSolver::solve()
{
	if(!isRunning()) //If a thread is already active, don't deal in unnecessary concurrency. Solution is fundamentally serial. 
		start(); //start() is inherited from QThread and executes run() in a new thread.
}

void ArmSolver::run()
{
	while(true)
	{
		if(solvesemaphore.available()<2) break; //We have to be able to snag two data points for interpolation in our ODE function.
		if(!solvesemaphore.tryAcquire(1,17)) break; //Try for at least 1 frame to start solving, in practice means single-threaded solution.
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
			es.clear();
			esdot.clear();
			etimes.clear();
			seeded=false;
			//Since we just emptied the queues, also empty the semaphore
			solvesemaphore.acquire(solvesemaphore.available());
			destructomutex.unlock();
			break; //We just emptied the semaphore
		}
		
		qst=point(y[0],y[1]); //Initial conditions for next run-through are final conditions from this one.
		qstdot=point(y[2],y[3]);
		
		if(model==TORQUESCALEDANDREFLEX)
		{
			point e, edot;
			etimes.push_back(tk);
			if(solveDes) {e=qm[1]-qst; edot=qmdot[1]-qstdot;}
			else {e=qst-qm[1]; edot=qstdot-qmdot[1];}
			es.push_back(e);
			esdot.push_back(edot);
			
			if (etimes.size()>2) //Only read elements that are actually there...
				while(etimes[1]<(tk-REFLEXDELAY)) //While the 0th point is no longer needed, ie. because we just used it
				{
					etimes.pop_front();
					es.pop_front();
					esdot.pop_front();
				}
		}
		
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
		
