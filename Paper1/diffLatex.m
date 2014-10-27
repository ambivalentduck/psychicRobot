function diffLatex(old,new)

%Add a new line after each period because of diff's fixation with lines.
system(['sed ''s/\. /\.\n/g'' ',old,' >o.txt']);
system(['sed ''s/\. /\.\n/g'' ',new,' >n.txt']);
system('diff -U 0 o.txt n.txt >diffed.txt'); %-U just makes this easier to parse. 0 means no context lines

oldcells=strsplit(fileread(old),'\n');
newstring=fileread(new);
changedlatex=newstring;
diffedstring=fileread('diffed.txt');
diffedstring=regexprep(diffedstring,'\@\@[+-\d,\s]*\@\@','%@*');
diffedcell=strsplit(diffedstring,'%@*');

red={};
blue={};

for k=1:length(diffedcell)
    stripped=strtrim(diffedcell{k});
    if length(stripped)<2
        continue
    end
    if (stripped(1)=='-' && stripped(2)~='-')||stripped(1)=='+'
        %Block of redactions then blocks of additions
        cellstrs=strsplit(diffedcell{k},'\n');
        removed={};
        added={};
        for kk=1:length(cellstrs)
            if isempty(cellstrs{kk})
                continue
            elseif cellstrs{kk}(1)=='-'
                removed{end+1}=regexprep(cellstrs{kk},'^-\s*','');
            else 
                added{end+1}=regexprep(cellstrs{kk},'^\+\s*','');
            end
        end
        red{k}=strjoin(removed(~isempty(removed)));
        blue{k}=strjoin(added(~isempty(added)));
    end
end

for k=1:length(diffedcell)
    %Four possibilities: 
    if isempty(red{k})&&isempty(blue{k}) %Only whitespace changed
        continue
    elseif ~isempty(red{k})&&~isempty(blue{k}) %Only substitutions
        changedlatex=strrep(changedlatex,blue{k},['\textcolor{red}{',red{k},'} \textcolor{blue}{',blue{k},'}']);
    elseif isempty(red{k}) %Only additions
        if length(blue{k})<3
            continue
        end
        changedlatex=strrep(changedlatex,blue{k},['\textcolor{blue}{',blue{k},'}']);
    else %Only lines removed.
        if isempty(strtrim(red{k}))
            continue
        end
        s=find(strcmp(oldcells,red{k}))
        if length(s)~=1
            continue
        end
        previous=oldcells{s-1}
        changedlatex=strrep(changedlatex,previous,[previous,' \textcolor{red}{',red{k},'}']);
    end
end
red'
blue'

fwrite(fopen('diffedlatex.tex','w'),changedlatex);
