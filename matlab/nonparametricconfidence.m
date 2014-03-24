clc
clear all

% Next is to align trajectories with 1.5mm bins
% --Throw away N*2.5% from the max and min extrema and note the min, max, mean for each bin.
% **This is a 95% confidence interval.
% --Plot this as a shaded region
% --Plot spaghetti that produced it on top
% --Plot median and mean on top of that.
%
% In theory, this can be done in inverse: what % of trials lie between a given point and the median?
% --Store lists by bin. So this is a list of lists.
% --For each point of interest, just ask what % of the appropriate list is less than that point.
% --Figure out how much significance it should take to be significant.

for k=2 %1:4
    load(['../Data/Data_pulse/pulse',num2str(k),'.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])

    clustermeS=zeros(length(trials)-1,1);
    clustermeE=clustermeS;
    for t=1:length(trials)-1
        plot(trials(t+1).x([1 end],1),trials(t+1).x([1 end],2),'.')
        clustermeS(t)=trials(t+1).x(1,1);
        clustermeE(t)=trials(t+1).x(end,1);
    end
    [cats,means]=kmeans([clustermeE; clustermeE],3,'emptyaction','singleton','start',[-0.175;-0.026;0.125]);
    means=sort(means);

    %For each direction/length/disturbance (2*2*5), plot an example from subject S
    %UNDER that example, plot the t-statistic: (Y-mean)/(sd(baseline)*sqrt(1+1/Nbaseline)
    % 95% confidence t<=.2

    nclean=2;
    clean=0*dcats;
    for kk=1:nclean
        clean=clean+[zeros(1,kk-1) dcats(1:end-(kk-1))];
    end
    clean=(clean==0)';

    figure(k)
    clf
    hold on

    for ST=1:3
        for EN=1:3
            if ST==EN
                continue
            end
            f=find((starts==ST)&(ends==EN)&clean);
            [ST EN length(f)]
            bins=[-inf min(means([ST EN])):.0015:max(means([ST EN])) inf];
            X=vertcat(trials(f).x);
            clear Tcat
            for kf=1:length(f)
                Tcat(kf).T=f(kf)+0*trials(f(kf)).x(:,1);
            end
            T=vertcat(Tcat.T);

            medy=zeros(length(bins)-1,1);
            lower=medy;
            upper=medy;

            for B=1:length(bins)-1
                inds=find((X(:,1)>=bins(B))&(X(:,1)<bins(B+1)));
                y=sort(X(inds,2));
                medy(B)=median(y);
                off=.025*length(y);
                if off<1
                    off=1;
                end
                frac=mod(off,1);
                yi=interp1(1:length(y),y,[off length(y)+1-off],'linear','extrap');
                lower(B)=yi(1);
                upper(B)=yi(2);
            end
            medy=smooth(medy);
            lower=smooth(lower);
            upper=smooth(upper);

            xoff=ST/3-means(ST);
            yoff=EN/25;

            H=bins(2:end-1);
            h=fill([H H(end:-1:1)]+xoff,[lower(2:end)' upper(end:-1:2)']+yoff,'k');
            set(h,'Facecolor',.5*[1 1 1]);

            nout=zeros(length(f),1);
            for kf=1:1:length(f)
                x=trials(f(kf)).x;
                %plot(x(:,1)+xoff,x(:,2)+yoff,'b.','markersize',.000001')
                for B=1:length(bins)-1
                    inds=find((x(:,1)>=bins(B))&(x(:,1)<bins(B+1)));
                    if sum(x(inds,2)<lower(B))||sum(x(inds,2)>upper(B))
                        nout(kf)=nout(kf)+1;
                    end
                end
            end
            nout2=sort(nout);
            red_lim=interp1(1:length(f),nout2,.95*length(f))

            for kf=1:1:length(f)
                x=trials(f(kf)).x;
                if nout(kf)>red_lim
                    color='r';
                else
                    color='b';
                    %continue
                end
                plot(x(:,1)+xoff,x(:,2)+yoff,[color,'.'],'markersize',.000001','linewidth',.000001')
            end
           
            plot(H+xoff,medy(2:end)+yoff,'k')
            plot(means(ST)+xoff,.5+yoff,'rx')
        end
    end
end
axis equal
axis off

