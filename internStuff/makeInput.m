clc
clear all

dist=.15;

center=[0 .48];
angles=(180+[30 150 270])'*pi/180;

targets=[center;center;center]+dist*[cos(angles) sin(angles)];

% #Each direction, force mag, EA gain
phases=[16 0 1;
    16 1.2 1;
    16 1.2 0;
    16 1.2 1;
    16 0 1];

q=0;
for P=1:size(phases,1)
    %Make a list of #1s, 2s and 3s
    deck=sort(mod(1:phases(P,1)*3,3)+1);
    deck=deck(randperm(phases(P,1)*3));
    for k=1:length(deck)
        q=q+1;
        REACHES(q).vals=[q targets(deck(k),1) targets(deck(k),2) 0 0 phases(P,2) -1 1 1 phases(P,3)];
        q=q+1;
        REACHES(q).vals=[q center(1) center(2) 0 0 phases(P,2) -1 1 1 phases(P,3)];

    end
end

reaches=vertcat(REACHES.vals)

fid=fopen('internInput.dat','w');
[stringy,ermsg]=sprintf('%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%g\n',reaches')
fprintf(fid,'%d\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%g\n',reaches')
fclose(fid);

%temptrial,&tempx,&tempy,&tempEarly,&tempLate,&tempWhite,&tempshape,&tempShowCursor, &tempAcquisitionsNeeded)