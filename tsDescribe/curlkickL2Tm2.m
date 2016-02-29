clc
clear all
close all

for SUB=1:6
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
    
    %isfirst=(order==2);
    %L2=L2(isfirst);
    %Tm2=Tm2(isfirst);
    
    f=find(~imag(L2)); %sanity check
    L2=L2(f);
    Tm2=Tm2(f);
    f=find(Tm2<50); %sanity check
    L2=L2(f);
    Tm2=Tm2(f);
    
    enstruct(SUB).Tm2=Tm2;
    [cumf,bins]=ecdf(Tm2);
    f=find(cumf<.8);
    [~,closest]=min(abs(cumf-.5));
    enstruct(SUB).Tm2F=cumf(f);
    enstruct(SUB).Tm2X=bins(f)/bins(closest);
    
    enstruct(SUB).L2oT2=L2.*Tm2;
    [cumf,bins]=ecdf(L2.*Tm2);
    f=find(cumf<.80);
    [~,closest]=min(abs(cumf-.5));
    enstruct(SUB).L2oT2F=cumf(f);
    enstruct(SUB).L2oT2X=bins(f)/bins(closest);
    
    nbins=15;
    sy=3;
    sx=2;
    
    continue 
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
    ylabel('log PDF(L^2T^{-2})')
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

figure(9)
clf
subplot(1,2,1)

X=vertcat(enstruct.L2oT2X);
X=[X ones(size(X))];
y=log(.85-vertcat(enstruct.L2oT2F));
%y=vertcat(enstruct.L2oT2F);

b=X\y;
yhat=X*b;

R2 = 1 - sum((y - yhat).^2)/sum((y - mean(y)).^2);
hold on
plot(X(:,1),y,'.','markersize',.75)
plot(X(:,1),X*b,'k-')
title(['R^2 = ',num2str(R2)])
ylabel('log PDF( L^2 T^{ -2} )')
xlabel('L^2 T^{ -2} Normalized to 50th percentile','interpreter','tex')
set(gca,'xtick',[0 1 2],'xticklabels',{'0','0.5','1'})
xlim([0 2])


subplot(1,2,2)
X=vertcat(enstruct.Tm2X);
X=[X ones(size(X))];
y=log(1.05-vertcat(enstruct.Tm2F));

b=X\y;
yhat=X*b;

R2 = 1 - sum((y - yhat).^2)/sum((y - mean(y)).^2);
hold on
plot(X(:,1),y,'.','markersize',.75)
plot(X(:,1),X*b,'k-')
title(['R^2 = ',num2str(R2)])
ylabel('log PDF( T^{ -2} )')
xlabel('T^{ -2} Normalized to 50th percentile','interpreter','tex')
set(gca,'yAxisLocation','right')
set(gca,'xtick',[0 1 2],'xticklabels',{'0','0.5','1'})
xlim([0 2])
