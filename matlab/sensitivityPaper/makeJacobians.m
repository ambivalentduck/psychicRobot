function makeJacobians

syms theta1 theta2 l1 l2 x0x x0y real
theta=[theta1; theta2];
x0=[x0x,x0y];
l=[l1; l2];

fk=fkinP(theta,l1,l2,x0);

fh=fopen('fJ.m','w');
fprintf(fh,'function out=fJ(theta,l1,l2)\ntheta1=theta(1);\ntheta2=theta(2);\n');
fprintf(fh,['fJ11=',vectorize(diff(fk(1),theta1)),';\n']);
fprintf(fh,['fJ21=',vectorize(diff(fk(2),theta1)),';\n']);
fprintf(fh,['fJ12=',vectorize(diff(fk(1),theta2)),';\n']);
fprintf(fh,['fJ22=',vectorize(diff(fk(2),theta2)),';\n']);
fprintf(fh,'out=[fJ11,fJ12;fJ21,fJ22];\n');
fclose(fh);

fh=fopen('getAlpha.m','w');
fprintf(fh,'function out=getAlpha(theta,omega,accel,l1,l2)\ntheta1=theta(1);\ntheta2=theta(2);\n');
fprintf(fh,['fJt1dx2=',vectorize(diff(fk(1),theta1,2)),';\n']);
fprintf(fh,['fJt1dy2=',vectorize(diff(fk(1),theta2,2)),';\n']);
fprintf(fh,['fJt2dx2=',vectorize(diff(fk(2),theta1,2)),';\n']);
fprintf(fh,['fJt2dy2=',vectorize(diff(fk(2),theta2,2)),';\n']);
fprintf(fh,['fJt1dxdy=',vectorize(diff(diff(fk(1),theta1),theta2)),';\n']);
fprintf(fh,['fJt2dxdy=',vectorize(diff(diff(fk(2),theta1),theta2)),';\n']);
fprintf(fh,'out=fJ(theta,l1,l2)\\(accel-[fJt1dx2*omega(1)+fJt1dxdy*omega(2),fJt1dy2*omega(2)+fJt1dxdy*omega(1);fJt2dx2*omega(1)+fJt2dxdy*omega(2),fJt2dy2*omega(2)+fJt2dxdy*omega(1)]*omega);\n');
fclose(fh);
