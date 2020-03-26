% --------------------------------------------------------------
% MAIN SCRIPT (2 public goods)
%
% This is the main script for running bacteria simulations 
% secreting two public goods.
% The code is used in the paper ''Evolution of specialization
% in dynamic fluids'' by G. Uppal and DC Vural
% Parameters are mostly chosen to have informative names,
% comments below should provide additional information. 
% To compare with values used in the paper, simulation units are
% scaled such that 1 unit of simulation time corresponds to 10^4 
% physical seeconds, and 1 unit of space corresponds to 1 mm. 
% Default values used below correspond to the (AND) values in 
% Table 1 of the paper.
% --------------------------------------------------------------


% === CLEAR WORKSPACE ===
close all; clear; clc;
% TIME SIMULATION:
tic; % start stopwatch
%----------------------------------------------------------------

% ===== USER DEFINED PARAMETERS =====

% === SPACE-TIME PARAMETERS ===
% SPACE:
xmin = -20;
xmax = 20;
ymin = -20;
ymax = 20;

% INITIALIZATION DOMAIN:
% A subdomain in which bacteria are randomly placed
% is given by the following parameters. By default,
% these are chosen to correspond to the whole simulation domain
xbmin = xmin;
xbmax = xmax;
ybmin = ymin;
ybmax = ymax;

% BOUNDARY CONDITIONS: (for chemicals and bacteria)
% (0 == periodic; 1 == no flux (Neumann); 
    % 2 == open -- option to have right end absorbing, not used in this paper )
xbnd = 0;
ybnd = 0;

% TIME:
run_time = 1.0;

% option to save intermediate times:
tsave = []; 

% === SYSTEM PARAMETERS ===
% diffusion:
db = 1.0;       % bacteria diffusion constant
d1 = 5;         % resource 1 diffusion constant
d2 = 5;         % resource 2 diffusion constant
dw = 15;        % waste diffusion constant

% == FITNESS PARAMETERS ==
% benefit/harm:
a1 = 0;        % resource benefit of just resource one -- not used in paper
a2 = 0; 

a12_or = 65;        % benefit constant in front of or term
a12_and = 0;        % benefit constant in front of and term 
                        % -- set one to zero to use just and or or
aw = 105;            % waste harm

% secretion rates:
s1 = 100;       % initial resource 1 secretion (subject to mutations)
s2 = 100;       % initial resource 2 secretion
sw = 90;       % initial waste secretion

% 0 == false, 1 == true:
s1rand = 0;     % option to start with random values for s1 
s2rand = 0;     % option to start with random values for s2
sbinary = 0;    % option to start with specialists where everyone 
                    % only secretes one of the two public goods

% secretion costs:
gamma1 = 0;     % constant cost for secretion
gamma2 = 0;     % constant cost for secretion -- not used

b1 = 0.1;       % resource 1 secretion cost
b2 = 0.1;       % resource 2 secretion cost

delta = 0.0;   % cost of secreting both

% saturation:
k1 = 0.01;      % resource 1 saturation constant
k2 = 0.01;      % resource 2 saturation constant
                    % these are for indiviually adding each good
                    % in a ``sum'' -- not used in paper
                    
kor = 0.01;     % resource saturation constant in OR fitness
k12 =  0.00003; % resource saturation constant in AND fitness
kw = 0.1;       % waste saturation constant

% chemical decays:
l1 = 50;        % resource 1 decay rate
l2 = 50;        % resource 2 decay rate
lw = 15;        % waste decay rate


% === FLOW VELOCITY PARAMETERS ===
% velocity_type options:
% 1 = Couette
% 2 = Hagen - Poisuille
% 3 = Rankine Vortex
% 4 = Constant (useful for sticky bacteria)
% else: no flow - (equivalent to setting vmax = 0)

velocity_type = 0;
vmax = 0;

% for vortex case:
rotation = 0;
vortex_radius = 0;


% === MUTATION PARAMETERS ===
mscale = 1e-7;

% mutation parameters:
mutProb1 =  2*mscale;   % mutation rate for first secretion function
mutDiff1 = 10;          % how much secretion rate changes due to a mutation
mdiffrand1 = 0;         % randomize strength of mutation between [-mdiff mdiff] 
                            % 0 == gaussian with variance mdiff
   % for second secretion function:
mutProb2 =  mutProb1; 
mutDiff2 = mutDiff1;    
mdiffrand2 = 0;

binaryMutation = 1;  % option = 1 means full mutation: either secrete or don't secrete

% Mutation scheme for specialization is to set binaryMutation = 1
% and set mutProb1 = mutProb2


% === BACTERIA INITIALIZATION ===
initialNumBacteria = 3500; 
initialNumGroups = 35;
aligngroups = 0;            % option to align groups along x axis 
                            % or to center a single group,
                            % otherwise groups are placed randomly in
                            % domain
                            
maxBacteria = 100000;       % maximum allowed bacteria
sticky = 0;                 % option to have only chemicals move with flow -- not used


% === SAVING AND VISUALIZATION ===
savePeriod = 10000;         % saves every this many time steps
saveChemicals = 0;          % option to save chemical data as well 
                                % rarely use -- will result in large output

% name of file for output:                                
savefile = sprintf('twogoodOR_T21_wmut_d1_%d_d2_%d_ds10_m20_g%.1f',...
    d1,d2,gamma1); 

graphics = 1;     % option to view simulation while running
framePeriod = 10; % display simulation every this many time steps
saveVideo = 0;    % option to save video (1 == true) -- graphics must be set to 1


%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
% ===== DERIVED AND RESERVED PARAMETERS =====
% (normally don't change these)

% === DISCRETIZATION ===
Nx = round((xmax-xmin)*5); % number of x nodes - 1, corresponding to dx = 0.2 
Ny = round((ymax-ymin)*5); % number of y nodes - 1, corresponding to dy = 0.2 
dt = 0.0001;                % time step

% === DELAYS ===
veldelay = 0.0020;       % delay time for flow 
mutdelay = 0.0020;       % delay for mutation 
repdelay = 0.0015;       % delay for reproduction

% === SAVE PARAMETERS ===
datechar = datestr(now,'mm-dd-yy_HH-MM-SS');
dataFile = sprintf('%s_%s_data.mat',savefile,datechar);
vidOption = 1;  % codec if saving video (2 == avi, else mp4)
vidFile = sprintf('%s_%s_vid',savefile,datechar);



%-------------------------------------------------------------------
%-------------------------------------------------------------------
%-------------------------------------------------------------------

% ===== RUN PROGRAM =====
[dataB, dataChem] = simulator(xmin,xmax,xbnd,Nx,ymin,ymax,ybnd,Ny,run_time,dt,...
    xbmin,xbmax,ybmin,ybmax, ...
    db,d1,d2,dw,a1,a2,a12_or,a12_and,aw,b1,b2,k1,k2,kor,k12,kw,s1,s2,sw,l1,l2,lw, ...
    gamma1,gamma2,delta, ...
    s1rand,s2rand,sbinary, ...
    velocity_type,vmax,rotation,vortex_radius,...
    mutProb1,mutDiff1,mdiffrand1,mutProb2,mutDiff2,mdiffrand2,binaryMutation,...
    sticky,initialNumBacteria,initialNumGroups,aligngroups,maxBacteria,...
    veldelay,mutdelay,repdelay,...
    savePeriod,saveChemicals,...
    graphics,framePeriod,saveVideo,vidOption,vidFile,...
    tsave,dataFile);

% === SAVE DATA ===
if saveChemicals == 1
    data = {dataB, dataChem};
else
    data = dataB;
end

parsave(dataFile,data);



% === DISPLAY RUN TIME ===
rntime = toc;
fprintf('\nProgram run time: %f \n',rntime);


%=======================================================================
% END OF FILE
%=======================================================================




