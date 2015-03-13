function start=onsetDetector(trial)

onset=find(vecmag(trial.v)>.1,1,'first');
start=max(onset-35,1);