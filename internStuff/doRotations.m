clc
clear all


subs=[2:9 23:30];

if ~exist('stiffness.mat','file')
    for S=subs
        if S<20
            prefix='intern';
        else
            prefix='output';
        end
        load(['./Data/',prefix,num2str(S),'.mat'])
        if ~isfield(trials,'xrot')||1
            coeffs=zeros(480,4);
            for T=2:480
%                             if (T>(2*96))&&(T<=(96*3))
%                                 p=trials(T).y;
%                             else
%                                 p=trials(T).x;
%                             end
                p=trials(T).x;
                x=[p(:,1)-trials(T).orig(1), p(:,2)-trials(T).orig(2)];
                trials(T).start=find(vecmag(x)>.02,1,'first');
                
                trials(T).xrot=rotateProgressError(trials(T).x,trials(T).orig,trials(T).targ);
                trials(T).frot=rotateProgressError(trials(T).f,trials(T).orig,trials(T).targ);
                trials(T).yrot=rotateProgressError(trials(T).y,trials(T).orig,trials(T).targ);
                reginds=trials(T).start+(0:39);
                mgt=mean(gradient(trials(T).t(reginds)));
                
                v_error=gradient(trials(T).xrot(reginds,2))/mgt;
                a_error=gradient(v_error)/mgt;
                %trials(T).regressX=[trials(T).xrot(reginds,2) v_error a_error]; %#ok<*SAGROW>
                trials(T).regressX=[trials(T).xrot(reginds,2) v_error a_error ones(size(a_error))]; %#ok<*SAGROW>
                
                fshift=-16; %80 ms
                if reginds(1)+fshift>=1
                    trials(T).regressF=trials(T).frot(reginds+fshift,2);
                else
                    undershoot=reginds(1)+fshift;
                    trials(T).regressF=[-trials(T-1).frot(end+undershoot:end,2);trials(T).frot(1:(reginds(end)+fshift),2)];
                end
                coeffs(T,:)=trials(T).regressX\trials(T).regressF;
            end
            save(['./Data/output',num2str(S),'.mat'],'trials','params')
        end
        inds=3:2:479;
        block=floor(inds/96)+1;
        coeffs=coeffs(inds,:);
        trls=trials(inds);
        ylabels={'Stiffness, N/m','Damping, N*s/m','Mass, kg'};
        b=zeros(4,5);
        
        %figure(S)
        %clf
        for B=1:5
            state=vertcat(trls(block==B).regressX);
            force=vertcat(trls(block==B).regressF);
            [b(:,B),bint{B}]=regress(force,state);
            %[b{B},bint{B}]=regress(force,state)
            %mdl=fitlm(state,force)
            for k=1:4
                %subplot(1,4,k)
                %hold on
                %plot(B*ones(sum(block==B),1)-.4*rand(sum(block==B),1),coeffs(block==B,k),'b.')
                if k==1
                    stiff(S,B)=median(coeffs(block==B,k));
                end
                %plot(B-.2,stiff(S,B),'rx')
            end
        end
    end
    save('stiffness.mat','stiff')
else
    load('stiffness.mat')
end

%% Post-hoc analysis

figure(50)
clf
subplot(2,4,1.5:3.5)
hold on
exp1=2:9;
exp2=23:30;
rn=.1*linspace(-1,1,8);

stiff=abs(stiff);

for s=1:8
    plot((1:5)-rn(s),stiff(exp1(s),:),'r-')
    plot((1:5)-rn(s),stiff(exp1(s),:),'r.')
    plot((1:5)-rn(s),stiff(exp2(s),:),'b-')
    plot((1:5)-rn(s),stiff(exp2(s),:),'b.')
end
set(gca,'xtick',1:5)
set(gca,'ytick',[0:50:200 225 250 275])
ylabs={};
for k=0:50:200
    ylabs{k/50+1}=num2str(k);
end
ylabs{end+1}='Cursor=Intent';
ylabs{end+1}='Cursor=Hand';
ylabs{end+1}='Forces On';
length(ylabs)
length([0:50:200 225 250 275])
set(gca,'yticklabel',ylabs)
ylim([0 300])

xlim([.5 5.5])
ylabel('Stiffness (N/m)')
xlabel('Block')
plot([1.5 4.5],275*[1 1],'k')
plot([2.5 3.5],225*[1 1],'k')
plot([.5 2.5],250*[1 1],'k',[3.5 5.5],250*[1 1],'k')

xl=[.7 2.3];
%Hypothesis 1: forces increase stiffness (2-1)
msize=10;
subplot(2,4,5)
miniplot(stiff(exp1,2),stiff(exp1,1),stiff(exp2,2),stiff(exp2,1),msize,xl,rn)
ylabel('\DeltaStiffness (N/m)')
title('2-1')

%Hypothesis 2: intent decreases stiffness despite forces (3-2)
subplot(2,4,6)
miniplot(stiff(exp1,3),stiff(exp1,2),stiff(exp2,3),stiff(exp2,2),msize,xl,rn)
title('3-2')

%Hypothesis 3: stiffness increases again when intent is removed (4-3)
subplot(2,4,7)
miniplot(stiff(exp1,4),stiff(exp1,3),stiff(exp2,4),stiff(exp2,3),msize,xl,rn)
title('4-3')

%Hypothesis 4: people are essentially adopting model stiffness
subplot(2,4,8)
load('modelstiff.mat')
ms=mean(modelstiff')';
miniplot(stiff(exp1,3),ms(exp1),stiff(exp2,3),ms(exp2),msize,xl,rn)
title('3-Model')