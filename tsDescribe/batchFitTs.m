clc
clear all

shift=zeros(8,1);
n=shift;
T=n;
costs=T;

for k=1:8
    [shift(k),n(k),T(k),costs(k)]=fitTs(k);
    %set(gcf,'position',[200   209   867   604]);
    %print('-dpng','-r300',['sub',num2str(k),'.png'])
end

[shift n T costs]

