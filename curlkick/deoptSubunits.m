function lumps_opt=deoptSubunits(lumps_init,t,y)

% F_VTR		"Value To Reach" (stop when ofunc < F_VTR)
F_VTR = 0.001; 

% I_D		number of parameters of the objective function 
I_D = length(lumps_init)/4; 
x_init=lumps_init(4*(1:I_D)); %Only the Ses

% FVr_minbound,FVr_maxbound   vector of lower and bounds of initial population
%    		the algorithm seems to work especially well if [FVr_minbound,FVr_maxbound] 
%    		covers the region where the global minimum is expected
%               *** note: these are no bound constraints!! ***
      FVr_minbound = x_init'*.5; 
      FVr_maxbound = x_init'*1.5; 
      I_bnd_constr = 0;  %1: use bounds as bound constraints, 0: no bound constraints      
            
% I_NP            number of population members
		I_NP = 15; 

% I_itermax       maximum number of iterations (generations)
		I_itermax = 1000; %Was 5000
       
% F_weight        DE-stepsize F_weight ex [0, 2]
		F_weight = 0.6; 

% F_CR            crossover probabililty constant ex [0, 1]
		F_CR = .6; 

% I_strategy     1 --> DE/rand/1:
%                      the classical version of DE.
%                2 --> DE/local-to-best/1:
%                      a version which has been used by quite a number
%                      of scientists. Attempts a balance between robustness
%                      and fast convergence.
%                3 --> DE/best/1 with jitter:
%                      taylored for small population sizes and fast convergence.
%                      Dimensionality should not be too high.
%                4 --> DE/rand/1 with per-vector-dither:
%                      Classical DE with dither to become even more robust.
%                5 --> DE/rand/1 with per-generation-dither:
%                      Classical DE with dither to become even more robust.
%                      Choosing F_weight = 0.3 is a good start here.
%                6 --> DE/rand/1 either-or-algorithm:
%                      Alternates between differential mutation and three-point-
%                      recombination.           

		I_strategy = 2;

% I_refresh     intermediate output will be produced after "I_refresh"
%               iterations. No intermediate output will be produced
%               if I_refresh is < 1
      I_refresh = 0;

% I_plotting    Will use plotting if set to 1. Will skip plotting otherwise.
      I_plotting = 0;

%-----Definition of tolerance scheme--------------------------------------
%-----The scheme is sampled at I_lentol points----------------------------
%-----tie all important values to a structure that can be passed along----
S_struct.t            = t;
S_struct.y            = y;
S_struct.lumps        = lumps_init;
S_struct.I_NP         = I_NP;
S_struct.F_weight     = F_weight;
S_struct.F_CR         = F_CR;
S_struct.I_D          = I_D;
S_struct.FVr_minbound = FVr_minbound;
S_struct.FVr_maxbound = FVr_maxbound;
S_struct.I_bnd_constr = I_bnd_constr;
S_struct.I_itermax    = I_itermax;
S_struct.F_VTR        = F_VTR;
S_struct.I_strategy   = I_strategy;
S_struct.I_refresh    = I_refresh;
S_struct.I_plotting   = I_plotting;

%********************************************************************
% Start of optimization
%********************************************************************
tic
[FVr_x,S_y,I_nf] = deopt('objSubunits',S_struct);
toc
final_cost=S_y

S=FVr_x;
lumps_opt=lumps_init;
lumps_opt(4*(1:I_D))=S;
lumps_opt=vec2lumps(lumps_opt);

