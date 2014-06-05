clc
clear all

figure(5)
clf
hold on

colors=[1 0 0; 0 1 0; 0 0 1; 1 1 0]; %RGBK

S=3;
ANECDOTES=ones(2*2*5,1);

vals=[];
for ST=1:3
    for EN=1:3
        if ST==EN
            continue
        end
        vals(end+1)=ST*10+EN;
    end
end
URC=unique(vals);

for k=1:4
    load(['../Data/Data_pulse/pulse',num2str(k),'W.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'Y.mat'])
    load(['../Data/Data_pulse/pulse',num2str(k),'I.mat'])

    %Start off by replicating the confidence interval code for disturbed movements.

    for ST=1:3
        for EN=1:3
            if ST==EN
                continue
            end
            direction=sign(EN-ST);
            distance=abs(EN-ST);
            RC=find((ST*10+EN)==URC);

            f=find((starts==ST)&(ends==EN)&dcats');
            bins=confidenceIntervals(ST,EN).bins;
            med=confidenceIntervals(ST,EN).med;
            upper=confidenceIntervals(ST,EN).upper;
            lower=confidenceIntervals(ST,EN).lower;

            outX=zeros(length(f),length(bins)-1);
            outY=zeros(length(f),length(bins)-1);

            for kf=1:length(f)
                %For each trial, just compare x and y via the bin find and record the boolean outcomes.
                X=trials(f(kf)).x;
                Y=trials(f(kf)).y;

                for B=1:length(bins)-1
                    indsX=find((X(:,1)>=bins(B))&(X(:,1)<bins(B+1)));
                    indsY=find((Y(:,1)>=bins(B))&(Y(:,1)<bins(B+1)));

                    if ~isempty(indsX)
                        [trash,i]=max(abs(X(indsX,2)-med(B)));
                        outX(kf,B)=(X(indsX(i),2)<lower(B))|(X(indsX(i),2)>upper(B));
                    else
                        outX(kf,B)=0;

                    end

                    if ~isempty(indsY)
                        [trash,i]=max(abs(Y(indsY,2)-med(B)));
                        outY(kf,B)=(Y(indsY(i),2)<lower(B))|(Y(indsY(i),2)>upper(B));
                    else
                        outY(kf,B)=0;
                    end
                end

                trials(f(kf)).outX=outX(kf,:);
                trials(f(kf)).outY=outY(kf,:);
            end
            %Now, by type of disturbance, foreach unique(dcat(f)), perform
            %kruskal wallis, which CAN have uneven group sizes.
            baseline=mean(confidenceIntervals(ST,EN).baseline,2);
            for D=1:5
                df=find(dcats(f)==D);
                disturbed=mean(outX(df,:),2);
                extracted=mean(outY(df,:),2);
                PKWhand(RC,D)=kruskalwallis([baseline; disturbed]',[ones(size(baseline)); 2*ones(size(disturbed))]','off');
                PKWextracted(RC,D)=kruskalwallis([baseline; extracted]',[ones(size(baseline)); 2*ones(size(extracted))]','off');
                if PKWextracted(RC,D)>.05
                    kruskalwallis([baseline; disturbed; extracted]',[ones(size(baseline)); 2*ones(size(disturbed)); 3*ones(size(extracted))]')
                end
            end



        end
    end
    PKWhand
    PKWextracted
end

%
%                     rc=find(urc==(10*starts(c)+ends(c)));
%                     xlookup=[0 .2 .6 1];
%                     xoff=(dir+3*dist)/3; %maps to 2 4 5 7
%                     yoff=(dcat-1)*.5;
%
%                     Y=trials(c).y(:,2);
%                     Y=Y-Y(1);
%                     X=trials(c).y(:,1);
%                     X=X-X(1);
%
%                     if (kk==1)&(k==1) %ANECDOTES(
%                         figure(5)
%                         plot(baselines(rc).edges(2:end-1)+xoff,baselines(rc).means(2:end)+yoff,'k')
%                         H=baselines(rc).edges(2:end-1);
%                         V1=baselines(rc).means(2:end)-1.96*baselines(rc).stds(2:end);
%                         V2=baselines(rc).means(2:end)+1.96*baselines(rc).stds(2:end);
%                         h=fill([H H(end:-1:1)]+xoff,[V1' V2(end:-1:1)']+yoff,'k');
%                         set(h,'Facecolor',.5*[1 1 1]);
%                         plot(X+xoff,Y+yoff,'r')
%                     end
%
%                     tstat=zeros(length(baselines(rc).edges)-1,1);
%
%                     for kkk=1:length(baselines(rc).edges)-1
%                         inds=find((X>=baselines(rc).edges(kkk))&(X<baselines(rc).edges(kkk+1)));
%                         if isempty(inds)
%                             tstat(kkk)=nan;
%                         else
%                             tstat(kkk)=abs(mean(Y(inds))-baselines(rc).means(kkk))/(baselines(rc).stds(kkk)); %*sqrt(1+1/baselines(rc).counts(kkk)));
%                         end
%                     end
%                     inds=find(isnan(tstat));
%                     for kkk=1:length(inds)
%                         iter=0;
%                         try
%                             while isnan(tstat(inds(kkk)-iter))
%                                 iter=iter+1;
%                             end
%                             lower=tstat(inds(kkk)-iter);
%                             Wl=1/iter;
%                         catch
%                             lower=0;
%                             Wl=0;
%                         end
%                         iter=0;
%                         try
%                             while isnan(tstat(inds(kkk)+iter))
%                                 iter=iter+1;
%                             end
%                             upper=tstat(inds(kkk)+iter);
%                             Wu=1/iter;
%                         catch
%                             upper=0;
%                             Wu=0;
%                         end
%                         tstat(inds(kkk))=(Wu*upper+Wl*lower)/(Wu+Wl);
%                     end
%
%                     figure(5)
%                     try
%                         forceinds=storeme(lookupStoreme(c,1),lookupStoreme(c,2)).forceson;
%                     catch
%                        forceinds=find(vecmag(trials(c).f)>.1);
%                     end
%                     %You have about 30cm spacing between anecdotes.
%                     subjecttop=-.3/4*(k-1)-.1;
%                     subjectheight=.3/4;
%                     rowtop=subjecttop-subjectheight*((kk-1)/length(F));
%                     rowbottom=subjecttop-subjectheight*((kk)/length(F));
%                     gap=.1*subjectheight/6;
%                     x=trials(c).x(:,1);
%                     x=x-x(1);
%                     h=fill([x(forceinds([1 end end 1]))]+xoff,[rowtop rowtop rowbottom rowbottom]+yoff,'k');
%                     set(h,'facecolor',[.7 1 .7]);
%
%                     edg=baselines(rc).edges;
%                     for kkk=2:length(edg)-2
%                         if tstat(kkk)>=3 %1.96
%                             h=fill(edg(kkk+[0 1 1 0])+xoff,[rowtop-gap/2 rowtop-gap/2 rowbottom+gap/2 rowbottom+gap/2]+yoff,'k');
%                             set(h,'facecolor',colors(k,:),'edgecolor',colors(k,:))
%                         end
%                     end
%                     if (sum(tstat)>0)&0
%                         outsidethelines=outsidethelines+1
%                         figure(678)
%                         xoff=0;
%                         yoff=outsidethelines/10;
%                         plot(baselines(rc).edges(2:end-1)+xoff,baselines(rc).means(2:end)+yoff)
%                         plot(baselines(rc).edges(2:end-1)+xoff,baselines(rc).means(2:end)-1.96*baselines(rc).stds(2:end)+yoff,'.')
%                         plot(baselines(rc).edges(2:end-1)+xoff,baselines(rc).means(2:end)+1.96*baselines(rc).stds(2:end)+yoff,'.')
%                         plot(X+xoff,Y+yoff,'r-o')
%                         plot(trials(c).x(:,1)-trials(c).x(1,1),trials(c).x(:,2)-trials(c).x(1,2)+yoff,'g-x')
%                         tst=tstat(2:end);
%                         tst(tst>.2)=.25;
%                         plot(baselines(rc).edges(2:end-1)+xoff,tst/10-.03+yoff,'k')
%                     end
%                 end
%
%             end
%         end
%     end
% end
% axis equal