function lumps=rulesFindLumps(t,y)

%Rule 1:
% Submovements touch end-to-end in time. This guarantees a stack depth of two or less, but not necessarily pure peaks.

%Rule 2:
% Submovements are perfect, straight-line 5th order polynomials.



