function arrow(x,y,color,p,varargin)

if nargin<5
    thick=1;
else
    thick=varargin{1};
end

headangle=15*(pi/180);

hold on
for k=1:size(x,1)
    %stem
    plot([x(k,1) y(k,1)],[x(k,2) y(k,2)],'color',color,'linewidth',thick)
    angle=atan2(y(k,2)-x(k,2),y(k,1)-x(k,1));
    mag=sqrt((y(k,2)-x(k,2)).^2+(y(k,1)-x(k,1)).^2);
    plot(y(k,1)+[p*mag*cos(pi-headangle+angle) 0 p*mag*cos(pi+headangle+angle)],y(k,2)+[p*mag*sin(pi-headangle+angle) 0 p*mag*sin(pi+headangle+angle)],'color',color,'linewidth',thick)
end

axis equal