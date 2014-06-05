clc
clear all

%something like...
c=0;

for k=1:4
    load('../Data/pulse',num2str(k),'.mat')

    f=find([trials.disturbcat]~=5);
    for kk=1:length(f)
        %Find onset of forces
        onset=0; %Find force magnitude > ? just plot it vs x on known pulse trials
        over=find(trials(f(kk)).t>=(trials(f(kk)).t(onset)+.15); %Sanity check to make sure forces are gone?
        c=c+1;
        %Have to convert forces and positions to torque and arm configurations
        %AND a reference arm configuration which can just have the y-component
        %zeroed out prior to conversion to joint angle space.

        [catme(c).predictedtorque catme(c).torque]=converteverything;
    end
end
