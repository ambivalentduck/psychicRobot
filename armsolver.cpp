#include "ode.h"

ArmSolver::ArmSolver(twoLinkArm::ArmParams P, bool solveQ)
{
	doQ=solveQ;
	arm=new twoLinkArm(P);
	voidpointer=(void*) this;
	sys = {statfunc, statjac, 4, voidpointer};
	driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,1e-6, 1e-6, 0.0);
}

ArmSolver::~ArmSolver()
{
	gsl_odeiv2_driver_free(driver);
}

point computeDynamics(point qDes, point qDesDot, point qDesDDot, point q, point qDot, point torque, mat2 Kp, mat2 Kd);

	

int ArmSolver::func(double t, const double y[], double f[])
{
	if(doQ)
	{
		point qDes=point(y[0],y[1]);
		point qDesDot=point(y[2],y[3]);
		f[0]=y[2];
		f[1]=y[3];
		point qDesDDot=computeInvDynamics(point q, point qDot, point qDDot, point qDes, point qDesDot, point torque, mat2 Kp, mat2 Kd);
	}
	return GSL_SUCCESS;
}

void ArmSolver::solve(double t1, double t2, point q1, point q1d, point q1dd, point q2,point q2d, point q2dd)
{
	
	
	
	
	
	if(N>0) {delete solved; delete times;}
	solved=new double[n];
	times=new double[n];
	N=n;
	solved[0]=y1;
	times[0]=t1;
	double t=t1;
	double y=y1;
	for(int k=1; k<=n; k++)
	{
		double tk=t1+(t2-t1)*double(k)/N;
		int status=gsl_odeiv2_driver_apply(driver, &t, tk, &y);
		times[k]=t;
		solved[k]=y;
		
		if(status!=GSL_SUCCESS) {std::cout<<"Oops. Cuh-rash at: "<<k<<std::endl; break;}
		
		std::cout<<times[k] TAB solved[k] << std::endl;
	}
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
		
		
		
