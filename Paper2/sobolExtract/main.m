clc
clear all

%% Just like with paper 1, we need to prepare a "ground" truth.

%Step 1: we need to presume some submovements and a production rule. Blindly.

ri=[0 .5];
rf=[.15 .5];
rt=.7;
t=[0:.005:2];

%If you presume you know the time onsets (ie 1 and 3 touch in the middle)
%and that the initial is position at time of onset
dcontrib=[.25 .5 .25];
angle=[-pi/6 0 pi/6];

(rf-ri)dcontrib(1)=cos(angle(1))*

xi=[ri; dcontrib(1)/2*
    
xf=[  ;  ; rf];

ol=.25;
ton=[0; ol*rt; (1-2*ol)*rt];
toff=[2*ol*rt; (1-ol)*rt; rt]; 
    
    