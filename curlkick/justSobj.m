function cost=justSobj(S)

global fit_t fit_y lumpsvec

copylumps=adjustLumpS(lumpsvec,S);

yhat=assembleLumps(fit_t',copylumps);

low_thresh=.075; %50 ms
high_thresh=1; %1 s
cost  = sum(sum((fit_y-yhat).^2))+100*sum((low_thresh-S).*(S<low_thresh))+10*sum((S-high_thresh).*(S>high_thresh));