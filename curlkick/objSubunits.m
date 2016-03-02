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
function S_MSE= objSubunits(FVr_temp, S_struct)
%---Check the passband first------------------------------------

t=S_struct.t;
y=S_struct.y;
lumps=S_struct.lumps;
lL=length(lumps)/4;
S=FVr_temp;
lumps(4*(1:lL))=S;

yhat=assembleLumps(t,lumps);

F_cost  = sum(sum((y-yhat).^2));

%   
%---End: tolerance scheme---------------------------------------
%----strategy to put everything into a cost function------------
S_MSE.I_nc      = lL;%S is positive
S_MSE.FVr_ca    = -S.*(S<=0)+(S-1).*(S>1); %Negative and low-frequency
S_MSE.I_no      = 1;%number of objectives (costs)
S_MSE.FVr_oa(1) = F_cost;