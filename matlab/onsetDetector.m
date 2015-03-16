function start=onsetDetector(trial)

onset=find(vecmag(trial.v)>.025,1,'first');
start=max(onset-35,1);