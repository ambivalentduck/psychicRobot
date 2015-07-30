function [vel,kerns,pos]=supMJ5P(w,tc,ts,t)

ldx=size(w,2);

vel=zeros(length(t),ldx);
kerns=zeros(length(t),length(tc));

if nargout>2
    pos=zeros(length(t),ldx);
end

for k=1:size(w,1)
    if nargout>2
        [kern,prog]=MJ5P(t,tc(k),ts(k));
    else
        kern=MJ5P(t,tc(k),ts(k));
    end
    kerns(:,k)=kern;
    
    vel=vel+kern*w(k,:);
    
    if nargout>2
        pos=pos+prog*w(k,:);
    end
end

%Test me: [~,~,blah]=supMJ5P([1 0;0 1;-1 0],[.25 .5 .75],.5*[1 1 1],0:.01:1);figure(1);clf;plot(blah(:,1),blah(:,2));