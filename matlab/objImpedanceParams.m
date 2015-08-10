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
function S_MSE= objImpedanceParams(FVr_temp, S_struct)

%---Check the passband first------------------------------------

qr = S_struct.qr;
qd = S_struct.qd;
tau_a = S_struct.tau_a;
l1 = S_struct.l1;

mpvec = FVr_temp(1:6);
kp0 = reshape(FVr_temp(7:10),2,2);
kp1 = reshape(FVr_temp(11:14),2,2);

tau_pred=zeros(size(tau_a));

for k=1:size(tau_a,1)
    [Dr,Cr]=computeDC_opt(qr(1:2,k),qr(3:4,k),l1,mpvec);
    [Dd,Cd]=computeDC_opt(qd(1:2,k),qd(3:4,k),l1,mpvec);

    Tff=Dd*qd(5:6,k)+Cd;
    Tin=Dr*qr(5:6,k)+Cr;
    Tmuscle=abs(Tin+tau_a(:,k));
    Ep=qr(1:2,k)-qd(1:2,k);
    Ev=qr(3:4,k)-qd(3:4,k);
    
    tau_pred(:,k)=Tff-Tin-(kp0+kp1*diag(Tmuscle))*(Ep+Ev/12);
end

F_cost  = sum(sum((tau_a-tau_pred).^2));

%   
%---End: tolerance scheme---------------------------------------
%----strategy to put everything into a cost function------------
S_MSE.I_nc      = 14;%no constraints
S_MSE.FVr_ca    = -FVr_temp.*(FVr_temp<0);%no constraint array
S_MSE.I_no      = 1;%number of objectives (costs)
S_MSE.FVr_oa(1) = F_cost;