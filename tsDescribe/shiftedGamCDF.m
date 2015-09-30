function F=shiftedGamCDF(X,U,N,A)
F=X-U;
F(F<=0)=0;
F(F>0)=gamcdf(F(F>0),N,A);