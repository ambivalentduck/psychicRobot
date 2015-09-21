clc
clear all
close all

for SUB=1:8
    load(['../Data/curlkick',num2str(SUB),'g.mat'])
    dt=mean(diff(trials(2).t));
    lt=length(trials);
    haslumps=zeros(lt,1);
    for k=1:lt
        haslumps(k)=~isempty(trials(k).nlumps);
    end
    f=find(haslumps) ;
    lumps=sum([trials.nlumps])
    order=zeros(lumps,1);
    L2=zeros(lumps,1);
    Tm2=zeros(lumps,1);
    L=0;
    for c=1:length(f)
        for cc=1:trials(f(c)).nlumps
            L=L+1;
            order(L)=cc;
            x=trials(f(c)).lumps(cc).y(:,1:2);
            L2(L)=sum((x(end,:)-x(1,:)).^2);
            Tm2(L)=(dt*size(x,1))^-2;
        end
    end
    
    %isfirst=isfirst==1;
    %L2=L2(isfirst);
    %Tm2=Tm2(isfirst);
    
    f=find(~imag(L2)); %sanity check
    L2=L2(f);
    Tm2=Tm2(f);
    f=find(Tm2<50); %sanity check
    L2=L2(f);
    Tm2=Tm2(f);
    
    nbins=15;
    sy=3;
    sx=2;
    
    figure(SUB)
    clf
    subplot(sy,sx,1)
    hold on
    [counts,bins]=hist(L2,nbins);
    plot(bins,log(counts),'.')
    [R2,A]=plotExp(bins,counts);
    title(['R^2 = ',num2str(R2),' A=',num2str(A)])
    ylabel('log PDF(L^2)')
    xlabel('L^2')
    subplot(sy,sx,2)
    ecdf(L2)
    ylabel('CDF(L^2)')
    xlabel('L^2')
    
    subplot(sy,sx,3)
    hold on
    [counts,bins]=hist(Tm2,nbins);
    plot(bins,log(counts),'.')
    [R2,A]=plotExp(bins,counts);
    title(['R^2 = ',num2str(R2),' A=',num2str(A)])
    ylabel('log PDF(T^{-2})')
    xlabel('T^{-2}')
    subplot(sy,sx,4)
    ecdf(Tm2)
    ylabel('CDF(T^{-2})')
    xlabel('T^{-2}')
    
    subplot(sy,sx,5)
    hold on
    [counts,bins]=hist(L2.*Tm2,nbins);
    plot(bins,log(counts),'.')
    [R2,A]=plotExp(bins,counts);
    title(['R^2 = ',num2str(R2),' A=',num2str(A)])
    ylabel('log PDF($\frac{L^2}{T^2}$)','interpreter','latex')
    xlabel('$\frac{L^2}{T^2}$','interpreter','latex')
    subplot(sy,sx,6)
    ecdf(L2.*Tm2)
    ylabel('CDF($\frac{L^2}{T^2}$)','interpreter','latex')
    xlabel('$\frac{L^2}{T^2}$','interpreter','latex')
    
    
    h=suplabel(['Curl-kick Experiment Submovement Decomposition, Subject ',num2str(SUB)],'t');
    set(h,'position',get(h,'position')-[0 .01 0 0])
    set(gcf,'position',[680   134   536   830])
    print('-dpng',['curlkick',num2str(SUB),'.png'])
    
end