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
function S_MSE= objMassStiff(FVr_temp, S_struct)

%---Check the passband first------------------------------------

Tmp = S_struct.Tmp;
Tsp = S_struct.Tsp;
Tnp = S_struct.Tnp;
FitInds = S_struct.FitInds;

alpha = FVr_temp(1);
beta = FVr_temp(2);
time_shift = 0; %round(FVr_temp(3));
time_shift=min(20,time_shift);
time_shift=max(-20,time_shift);

Tnp_sim = alpha*Tmp(FitInds,:) + beta*Tsp(FitInds,:);

F_cost  = sum(sum((Tnp_sim - Tnp(FitInds+time_shift,:)).^2)); % + sum((Tnp_sim(:,2) - Tnp(floor(time_index)+1:end,2)).^2);

%   
%---End: tolerance scheme---------------------------------------
%----strategy to put everything into a cost function------------
S_MSE.I_nc      = 0;%no constraints
S_MSE.FVr_ca    = 0;%no constraint array
S_MSE.I_no      = 1;%number of objectives (costs)
S_MSE.FVr_oa(1) = F_cost;