clc
clear all

n=5:8;

% workflow={'','addSubjectPulse';
%     'W','findKpGain';
%     'Y','extractPulses';
%     'I','nonparametricconfidence';
%     'IT','nonparametricconfidenceTime';
%     'NEVERMATCH','finalfig4'
%     };

workflow={'','addSubjectPulse';
    'W','findKpGain';
    'Y','extractPulses';
    'U','extractUndisturbed';
    };



mfiles=dir('.');
datafiles=dir('../Data/Data_pulse/');

for k=n
    flag=0;
    for kk=1:size(workflow,1)
        if (~exist(['../Data/Data_pulse/pulse',num2str(k),workflow{kk,1},'.mat'],'file'))||flag
            feval(workflow{kk,2},k);
            flag=1;
        elseif mfiles(strcmp({mfiles.name},[workflow{kk,2},'.m'])).datenum > datafiles(strcmp({datafiles.name},['pulse',num2str(k),workflow{kk,1},'.mat'])).datenum
            feval(workflow{kk,2},k);
            flag=1;
        end
    end
end