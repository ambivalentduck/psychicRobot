function start=curlKickOnsetDetector(trial)

onset=find(vecmag(trial.v)>.1,1,'first');
start=max(onset-10,1);

%[~,peaks]=findpeaks(vecmag(trial.v),'MinpeakHeight',.3);
%start=max(1,peaks(1)-40);