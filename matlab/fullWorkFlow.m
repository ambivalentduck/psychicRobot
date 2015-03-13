clc
clear all

n=1:8;

%Add differential optimization just in case it's needed
addpath ../../DeOpt/

% workflow={'','addSubjectPulse';
%     'W','findKpGain';
%     'Y','extractPulses';
%     'I','nonparametricconfidence';
%     'IT','nonparametricconfidenceTime';
%     'NEVERMATCH','finalfig4'
%     };

% workflow={'','addSubjectPulse';
%     'W','findWhiteKpGain';
%     'Y','extractPulses'
%     };

workflow={'W','findWhiteKpGain';
    'Y','extractPulses'
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