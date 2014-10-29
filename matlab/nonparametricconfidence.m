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

NBINSOUT=3;

for k=4 %1:4
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'])

    clustermeS=zeros(length(trials)-1,1);
    clustermeE=clustermeS;
    for t=1:length(trials)-1
        clustermeS(t)=trials(t+1).x(1,1);
        clustermeE(t)=trials(t+1).x(end,1);
    end
    [cats,means]=kmeans([clustermeE; clustermeE],3,'emptyaction','singleton','start',[-0.175;-0.026;0.125]);
    means=sort(means);

    nclean=2;
    dcats=[trials.disturbcat];
    starts=[trialInfo.startcat];
    ends=[trialInfo.endcat];
    
    clean=0*dcats;
    for kk=1:nclean
        clean=clean+[zeros(1,kk-1) dcats(1:end-(kk-1))];
    end
    clean=(clean==0);

    figure(k)
    clf
    hold on


    figure(k+4)
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
                Tinds=T(inds);
                Yinds=X(inds,2);
                medYinds=median(Yinds);
                uT=unique(Tinds);
                y=0*uT;
                for uTk=1:length(uT)
                    finds=find(Tinds==uT(uTk));
                    [trash,i]=max(abs(Yinds(finds)-medYinds));
                    y(uTk)=Yinds(finds(i));
                end
                y=sort(y);
                medy(B)=median(y);
                off=.025*length(y);
                if off<2
                    off=2;
                end
                yi=interp1(1:length(y),y,[off length(y)+1-off],'linear','extrap');
                lower(B)=yi(1);
                upper(B)=yi(2);
            end
            medy=smooth(medy);
            lower=smooth(lower);
            upper=smooth(upper);

            xoff=ST/3-means(ST);
            yoff=EN/25;

            figure(k)
            H=bins(2:end-1);
            h=fill([H H(end:-1:1)]+xoff,[lower(2:end)' upper(end:-1:2)']+yoff,'k');
            set(h,'Facecolor',.5*[1 1 1]);

            nout=zeros(length(f),length(bins)-1);
            for kf=1:1:length(f)
                x=trials(f(kf)).x;
                for B=1:length(bins)-1
                    inds=find((x(:,1)>=bins(B))&(x(:,1)<bins(B+1)));
                    nout(kf,B)=sum(x(inds,2)<lower(B))||sum(x(inds,2)>upper(B));
                end
            end
            nout2=sum(nout,2);
            red_lim=interp1(1:length(f),sort(nout2),.95*length(f));

            for kf=1:1:length(f)
                x=trials(f(kf)).x;
                if nout2(kf)>red_lim
                    color='r';
                else
                    color='b';
                end
                plot(x(:,1)+xoff,x(:,2)+yoff,[color,'.'],'markersize',.000001','linewidth',.000001')
            end

            confidenceIntervals(ST,EN).med=medy;
            confidenceIntervals(ST,EN).lower=lower;
            confidenceIntervals(ST,EN).upper=upper;
            confidenceIntervals(ST,EN).bins=bins;
            confidenceIntervals(ST,EN).baseline=nout;

            plot(H+xoff,medy(2:end)+yoff,'k')
            plot(means(ST)+xoff,.5+yoff,'gx')
            axis equal
            axis off

            figure(k+4)
            for IT=1:5
                H=bins(2:end-1);
                h=fill([H H(end:-1:1)]+xoff,[lower(2:end)' upper(end:-1:2)']+yoff+(IT-1)/4,'k');
                set(h,'Facecolor',.5*[1 1 1]);
            end

            f=find((starts==ST)&(ends==EN)&dcats);
            for kf=1:length(f)
                %for each trial need 2x2 things: N-point moving nout, sum
                %nout, BUT both start with bool in/outside the boundaries.
                x=trials(f(kf)).x;
                y=trials(f(kf)).y;
                if isempty(y)
                        continue
                    end
                trials(f(kf)).xoutB=zeros(length(bins)-1,1);
                trials(f(kf)).youtB=zeros(length(bins)-1,1);

                for B=1:length(bins)-1
                    inds=find((x(:,1)>=bins(B))&(x(:,1)<bins(B+1)));
                    if isempty(inds)&&(B>1)&&(B<length(bins)-1)
                        %When a bin is missing, just snag the closest point.
                        [trash,inds]=min(abs(x-(bins(B)+bins(B+1))/2));
                    end
                    trials(f(kf)).xoutB(B)=sum(x(inds,2)<lower(B))||sum(x(inds,2)>upper(B));

                    inds=find((y(:,1)>=bins(B))&(y(:,1)<bins(B+1)));
                    if isempty(inds)&&(B>1)&&(B<length(bins)-1)
                        %When a bin is missing, just snag the closest point.
                        [trash,inds]=min(abs(y-(bins(B)+bins(B+1))/2));
                    end
                    trials(f(kf)).youtB(B)=sum(y(inds,2)<lower(B))||sum(y(inds,2)>upper(B));
                end

                trials(f(kf)).xout=zeros(size(x,1),1);
                for kx=1:size(x,1)
                    B=find((x(kx,1)>=bins(1:end-1))&(x(kx,1)<bins(2:end)),1,'first');
                    trials(f(kf)).xout(kx)=(x(kx,2)<lower(B))||(x(kx,2)>upper(B));
                end
                firstxout=find(cumsum(trials(f(kf)).xout)>=NBINSOUT,1,'first');

                trials(f(kf)).yout=zeros(size(y,1),1);
                for ky=1:size(y,1)
                    B=find((y(ky,1)>=bins(1:end-1))&(y(ky,1)<bins(2:end)),1,'first');
                    trials(f(kf)).yout(ky)=(y(ky,2)<lower(B))||(y(ky,2)>upper(B));
                end
                firstyout=find(cumsum(trials(f(kf)).yout)>=NBINSOUT,1,'first');

                yoff2=yoff+(dcats(f(kf))-1)/4
                
                plot(x(:,1)+xoff,x(:,2)+yoff2,'k')
                qscale=1000;
                dsample=3;
                quiver(x(1:dsample:end,1)+xoff,x(1:dsample:end,2)+yoff2,trials(f(kf)).f(1:dsample:end,1)/qscale,trials(f(kf)).f(1:dsample:end,2)/qscale,0,'g')
                plot(x(firstxout:end,1)+xoff,x(firstxout:end,2)+yoff2,'k.')

                plot(y(:,1)+xoff,y(:,2)+yoff2,'r')
                plot(y(firstyout:end,1)+xoff,y(firstyout:end,2)+yoff2,'r.')
                axis equal
                axis off
            end
        end
    end

    %Add the analysis from finalfig4 here. Or maybe just above. Plot trials
    %and show when a 3-5 frame moving average puts them out maybe by
    %switching from . to -


    save(['../Data/Data_pulse/pulse',num2str(k),'I.mat'],'confidenceIntervals')
end