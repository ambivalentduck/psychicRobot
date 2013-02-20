#include "ode.h"

ArmSolver::ArmSolver(twoLinkArm::ArmParams P, double tspacing, bool solveIntent, bool constImpedance)
{
	solveDes=solveIntent;
	tSpacing=tspacing;
	arm=new twoLinkArm(P);
	voidpointer=(void*) this;
	sys = {statfunc, statjac, 4, voidpointer};
	driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,1e-6, 1e-6, 0.0);
}

ArmSolver::~ArmSolver()
{
	gsl_odeiv2_driver_free(driver);
	delete arm;
}

int ArmSolver::func(double t, const double y[], double f[])
{
	double s=(t-t[0])/(t[1]-t[0]);
	double oms=1l-s;
	point qmi=qm[0]*oms+qm[1]*s;
	point qmdoti=qmdot[0]*oms+qmdot[1]*s;
	point qmddoti=qmddot[0]*oms+qmddot[1]*s;
	point torquemi=torquem[0]*oms+torquem[1]*s;
	
	mat2 kpi;
	mat2 kdi;
	if(constImpedance) {kpi=Kp; kdi=Kd;}
	else {kpi=Kpm[0]*oms+Kpm[1]*s; kdi=Kdm[0]*oms+Kdm[1]*s;}
	
	f[0]=y[2];
	f[1]=y[3];
	
	point qsDDot;
	
	if(solveDes) qsDDot=arm->computeInvDynamics(qmi, qmdoti, qmddoti, point(y[0],y[1]), point(y[2],y[3]), torquemi, kpi, kdi);
	else qsDDot=arm->computeDynamics(qmi, qmdoti, qmddoti, point(y[0],y[1]), point(y[2],y[3]), torquemi, kpi, kdi);
	
	f[2]=qsDDot[0];
	f[3]=qsDDot[1];
	
	return GSL_SUCCESS;
}

void ArmSolver::push(double t, point p, point v, point a)
{
	point q;
	while(arm->ikin(p,q)) arm->moveShoulder(point(0,-.01));
	{
		
		
	point q = ikin(p);
	
	
	
	
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
		
		
		
