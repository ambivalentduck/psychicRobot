clc
clear all
close all

outs=zeros(10,4);

mcount=0;
for k=[1:2 4:20 22]
    subname=num2str(k,'%2.2d');
    fname=['./matData/3Dfree_exp_',subname,'stroke.mat'];

    if ~exist(fname,'file')||1
        %outputsS06D3F3_org_session
        load(['../org_data/outputsS',subname,'D1F1_org_session.mat']);
        dt=.01; %mean(gradient(vertcat(Pointer_arm.time)));
        x=vertcat(Pointer_arm.smooth_pos);
        v=[gradient(x(:,1)) gradient(x(:,2)) gradient(x(:,3))]/dt;
        a=dt*[gradient(v(:,1)) gradient(v(:,2)) gradient(v(:,3))]/dt;
        t=0:dt:dt*(size(x,1)-1);
        timedomain.pre.pointer.x=x;
        timedomain.pre.pointer.v=v;
        timedomain.pre.pointer.a=a;
        timedomain.pre.pointer.t=t';

        dt=mean(gradient(vertcat(Robot_arm.time)));
        x=vertcat(Robot_arm.smooth_pos);
        v=[gradient(x(:,1)) gradient(x(:,2)) gradient(x(:,3))]/dt;
        a=dt*[gradient(v(:,1)) gradient(v(:,2)) gradient(v(:,3))]/dt;
        timedomain.pre.robot.x=x;
        timedomain.pre.robot.v=v;
        timedomain.pre.robot.a=a;
        timedomain.pre.robot.t=t';

        load(['../org_data/outputsS',subname,'D6F3_org_session.mat']);
        dt=.01; %mean(gradient(vertcat(Pointer_arm.time)));
        x=vertcat(Pointer_arm.smooth_pos);
        v=[gradient(x(:,1)) gradient(x(:,2)) gradient(x(:,3))]/dt;
        a=dt*[gradient(v(:,1)) gradient(v(:,2)) gradient(v(:,3))]/dt;
        t=0:dt:dt*(size(x,1)-1);
        timedomain.post.pointer.x=x;
        timedomain.post.pointer.v=v;
        timedomain.post.pointer.a=a;
        timedomain.post.pointer.t=t';

        dt=mean(gradient(vertcat(Robot_arm.time)));
        x=vertcat(Robot_arm.smooth_pos);
        v=[gradient(x(:,1)) gradient(x(:,2)) gradient(x(:,3))]/dt;
        a=dt*[gradient(v(:,1)) gradient(v(:,2)) gradient(v(:,3))]/dt;
        timedomain.post.robot.x=x;
        timedomain.post.robot.v=v;
        timedomain.post.robot.a=a;
        timedomain.post.robot.t=t';

        save(fname,'timedomain');
    else
        load(fname)
    end

    figure(1)
    clf

    for p={'pre','post'}
        for pp={'robot','pointer'}
            clf
            m=boltzmannsubjectplots(timedomain.(p{1}).(pp{1}).t,timedomain.(p{1}).(pp{1}).x,timedomain.(p{1}).(pp{1}).v,timedomain.(p{1}).(pp{1}).a);
            mcount=mcount+1;
            metrics(mcount).mets=m;
            metrics(mcount).prepost=p{1};
            metrics(mcount).hand=pp{1};
            metrics(mcount).sub=k;
            
            suplabel(['Stroke Subject ',subname,p{1},' ',pp{1}],'t')
            set(gcf,'position',[76 11 1195 925])
            print('-dtiff','-r300',['./figs/3Dsummary',subname,p{1},'_',pp{1},'_stroke.tiff'])
            print('-dpng','-r300',['./figs/3Dsummary',subname,p{1},'_',pp{1},'_stroke.png'])
        end
    end
end

save('3dmetrics.mat','metrics')