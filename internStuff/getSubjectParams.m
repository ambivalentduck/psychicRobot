function params=getSubjectParams(num)

%Every other case is a new subject
switch num
    case {'2'} 
        l1=.36;
        l2=.38;
        mass=145;
        shoulder=[.06 .95];
    case {'23'} 
        l1=.36;
        l2=.38;
        mass=145;
        shoulder=[.03 .97];
        
    case {'3'} 
        l1=.35;
        l2=.36;
        mass=145;
        shoulder=[.04 .935];
    case {'24'} 
        l1=.35;
        l2=.36;
        mass=145;
        shoulder=[0 .93];
        
    case {'4'} 
        l1=.33;
        l2=.38;
        mass=180;
        shoulder=[.038 .93];
    case {'25'} 
        l1=.33;
        l2=.38;
        mass=180;
        shoulder=[0 .95];
        
    case {'5'} 
        l1=.34;
        l2=.38;
        mass=187;
        shoulder=[.0635 .97];
    case {'26'} 
        l1=.34;
        l2=.38;
        mass=187;
        shoulder=[.01 .97];  
    
    %Only case 6 for this subject    
    case {'6'} 
        l1=.317;
        l2=.3048;
        mass=135;
        shoulder=[.0762 .95];
        
    case {'7'} 
        l1=.3175;
        l2=.3556;
        mass=230;
        shoulder=[.0762 .988];
    case {'27'} 
        l1=.3175;
        l2=.3556;
        mass=230;
        shoulder=[0 .97];   
        
    %Only case 8 for this subject
    case {'8'} 
        l1=.381;
        l2=.3937;
        mass=195;
        shoulder=[.08255 1.06];
    %Only case 8 for this subject
    case {'9'} 
        l1=.33;
        l2=.3175;
        mass=128;
        shoulder=[.0508 1];
        
        
    case {'christine','christine2','christine21'}
        l1=.33;
        l2=.34;
        mass=120;
        shoulder=[.1 .96]; 
    case {'tejas','output22'}
        %these are for the output22 data!!!
        l1= .36;
        l2= .33;
        mass = 140;
        shoulder=[0 .92];
end

params.l1=l1;
params.l2=l2;
params.shoulder=shoulder;

if exist('mass','var')
    params.mass=0.453592*mass;
else
    params.mass=((l1+l2)/.77)*91; %Shitty locally linear scaling of weight with height, remember that 200 lbs = 91 kg.
end