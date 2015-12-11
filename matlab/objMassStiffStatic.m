%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:         S_MSE= objfun(FVr_temp, S_struct)
% Author:           Rainer Storn
% Description:      Implements the cost function to be minimized.
% Parameters:       FVr_temp     (I)    Paramter vector
%                   S_Struct     (I)    Contains a variety of parameters.
%                                       For details see Rundeopt.m
% Return value:     S_MSE.I_nc   (O)    Number of constraints
%                   S_MSE.FVr_ca (O)    Constraint values. 0 means the constraints
%                                       are met. Values > 0 measure the distance
%                                       to a particular constraint.
%                   S_MSE.I_no   (O)    Number of objectives.
%                   S_MSE.FVr_oa (O)    Objective function values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S_MSE= objMassStiffStatic(FVr_temp, S_struct)
%---Check the passband first------------------------------------

X=S_struct.X;
Y=S_struct.Y;

c = FVr_temp(1);
L1 = FVr_temp(2);
L2 = FVr_temp(3);
x0 = FVr_temp(4:5);
mass = FVr_temp(6);

XfromS=vecmag([X(:,1)-x0(1) X(:,2)-x0(2)]);
YfromS=vecmag([Y(:,1)-x0(1) Y(:,2)-x0(2)]);
lArm=L1+L2;
minArm=abs(L1-L2);
if any((XfromS>lArm)|(YfromS>lArm)|(XfromS<minArm)|(YfromS>minArm))
    F_cost=inf;
else
    [TBMP,TSP,TNP]=cart2modelStatic(X,Y,L1,L2,x0,mass);
    Tnp_sim = TBMP + c*TSP;
    F_cost  = sum(sum((Tnp_sim - TNP).^2));
end

%   
%---End: tolerance scheme---------------------------------------
%----strategy to put everything into a cost function------------
S_MSE.I_nc      = 0;%no constraints
S_MSE.FVr_ca    = 0;%no constraint array
S_MSE.I_no      = 1;%number of objectives (costs)
S_MSE.FVr_oa(1) = F_cost;