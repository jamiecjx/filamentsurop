%-------------------------------------------------------------------------
%  Package up data in a .mat file depending on parameters, for use in JFNK
%-------------------------------------------------------------------------
global new_x	% Current best x
global epsJ		% epsilon used in Jacobian approximation
global ndts		% Number of timesteps taken in period T
global fixT		% Fix T for equilibrium, rather than PO solution
global f		% Follower force (parameter of dynamical system)
global N        % Number of segments
global FFTip    % True if simulation is for a follower force at the tip
global FFLength % True if simulation is for a follower force along the filament length
global Nf       % Number of filaments
global d        % Distance between filaments

% Simulation parameters (to edit)
Nf = 2;
N = 20;
d = 2.2*N;
ndts = 200;
t0 = 9.65;
x0 = u2(end,:);
new_x = [t0, x0]';
f = 100;

FFTip = true;
FFLength = false;

n       = 3*(N-1)*Nf+1;	% Dimension of system, including unknown params
mgmres  = 5;	% max GMRES iterations
nits    = 150 ;	    % max Newton iterations
rel_err = 1d-6;%1d-8;	% Relative error |F|/|x| %% i've changed this to be the error |F|: edited bottom of this script, and saveorbit

del     = -1d0 ;	% These rarely need changing for any problem
mndl    = 1d-20 ;
mxdl    = 1d+20 ;
gtol    = 1d-4;%1d-3 ;
epsJ    = 1d-6;%1d-5 ;

fixT = 0;

save("JFNKinput.mat", "Nf", "N", "d", "ndts", "t0", "x0", "f", ...
    "FFTip", "FFLength", "mgmres", "nits", "rel_err", "del", "mndl", ...
    "mxdl", "gtol", "epsJ", "fixT")