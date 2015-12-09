function start=curlKickOnsetDetector(trial)

onset=find(vecmag(trial.v)>.1,1,'first');
start=max(onset-10,1);