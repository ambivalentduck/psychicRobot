clc
clear all

SUB=3;
figure(1)
clf
hold on

rs=processPulse(SUB);
for u=unique([rs.reachcat])
    f=find([rs.reachcat]==u);
    figure(u+10)
    clf
    hold on
    for k=1:length(f)
        rs(f(k)).y=rotateProgressError(rs(f(k)).x,rs(f(k)).x0,rs(f(k)).xf);
        plot3(rs(f(k)).y(:,1),rs(f(k)).y(:,2),vecmag(rs(f(k)).v).^2,'.')
    end
    y=vertcat(rs(f).y);
    yp=hist3(y,[20 20]);
    %surf(log(yp))
    
    figure(1)
end
        