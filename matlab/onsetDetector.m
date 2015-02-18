function start=onsetDetector(trial)

onset=find(vecmag(trial.v)>.05,1,'first');
start=max(onset-35,1);