function [fJinline, getAlphainline, getAccelinline]=makeJacobians(varargin)

if nargin>0
    alphaName=varargin{1};
    fJName=varargin{2};
end

syms theta1 theta2 real
theta=[theta1; theta2];

fk=fkin(theta);

fJ11=inline(vectorize(diff(fk(1),theta1)));
fJ21=inline(vectorize(diff(fk(2),theta1)));
fJ12=inline(vectorize(diff(fk(1),theta2)));
fJ22=inline(vectorize(diff(fk(2),theta2)));

fJt1dx2=inline(vectorize(diff(fk(1),theta1,2)));
fJt1dy2=inline(vectorize(diff(fk(1),theta2,2)));
fJt2dx2=inline(vectorize(diff(fk(2),theta1,2)));
fJt2dy2=inline(vectorize(diff(fk(2),theta2,2)));
fJt1dxdy=inline(vectorize(diff(diff(fk(1),theta1),theta2)));
fJt2dxdy=inline(vectorize(diff(diff(fk(2),theta1),theta2)));

getAccelinline=@(theta,omega,alpha)[fJ11(theta(1),theta(2)),fJ12(theta(1),theta(2));fJ21(theta(1),theta(2)),fJ22(theta(1),theta(2))]*alpha+[fJt1dx2(theta(1),theta(2))*omega(1)+fJt1dxdy(theta(1),theta(2))*omega(2),fJt1dy2(theta(1),theta(2))*omega(2)+fJt1dxdy(theta(1),theta(2))*omega(1);fJt2dx2(theta(1),theta(2))*omega(1)+fJt2dxdy(theta(1),theta(2))*omega(2),fJt2dy2(theta(1),theta(2))*omega(2)+fJt2dxdy(theta(1),theta(2))*omega(1)]*omega;

fJinline=@(theta) [fJ11(theta(1),theta(2)),fJ12(theta(1),theta(2));fJ21(theta(1),theta(2)),fJ22(theta(1),theta(2))];

if nargin>0
    fh=fopen([fJName,'.m'],'w');
    fprintf(fh,['function out=',fJName,'(theta)\ntheta1=theta(1);\ntheta2=theta(2);\n']);
    fprintf(fh,['fJ11=',vectorize(diff(fk(1),theta1)),';\n']);
    fprintf(fh,['fJ21=',vectorize(diff(fk(2),theta1)),';\n']);
    fprintf(fh,['fJ12=',vectorize(diff(fk(1),theta2)),';\n']);
    fprintf(fh,['fJ22=',vectorize(diff(fk(2),theta2)),';\n']);
    fprintf(fh,'out=[fJ11,fJ12;fJ21,fJ22];\n');
    fclose(fh);
end

getAlphainline=@(theta,omega,accel) fJ(theta)\(accel-[fJt1dx2(theta(1),theta(2))*omega(1)+fJt1dxdy(theta(1),theta(2))*omega(2),fJt1dy2(theta(1),theta(2))*omega(2)+fJt1dxdy(theta(1),theta(2))*omega(1);fJt2dx2(theta(1),theta(2))*omega(1)+fJt2dxdy(theta(1),theta(2))*omega(2),fJt2dy2(theta(1),theta(2))*omega(2)+fJt2dxdy(theta(1),theta(2))*omega(1)]*omega);
if nargin>0
    fh=fopen([alphaName,'.m'],'w');
    fprintf(fh,['function out=',alphaName,'(theta,omega,accel)\ntheta1=theta(1);\ntheta2=theta(2);\n']);
    fprintf(fh,'global fJ\n');
    fprintf(fh,['fJt1dx2=',vectorize(diff(fk(1),theta1,2)),';\n']);
    fprintf(fh,['fJt1dy2=',vectorize(diff(fk(1),theta2,2)),';\n']);
    fprintf(fh,['fJt2dx2=',vectorize(diff(fk(2),theta1,2)),';\n']);
    fprintf(fh,['fJt2dy2=',vectorize(diff(fk(2),theta2,2)),';\n']);
    fprintf(fh,['fJt1dxdy=',vectorize(diff(diff(fk(1),theta1),theta2)),';\n']);
    fprintf(fh,['fJt2dxdy=',vectorize(diff(diff(fk(2),theta1),theta2)),';\n']);
    fprintf(fh,'out=fJ(theta)\\(accel-[fJt1dx2*omega(1)+fJt1dxdy*omega(2),fJt1dy2*omega(2)+fJt1dxdy*omega(1);fJt2dx2*omega(1)+fJt2dxdy*omega(2),fJt2dy2*omega(2)+fJt2dxdy*omega(1)]*omega);\n');
    fclose(fh);
end