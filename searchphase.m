function u = searchphase(fa, da)

global new_x	% Current best x
global ndts		% Number of timesteps taken in period T
global fixT		% Fix T for equilibrium, rather than PO solution
global f		% Follower force (parameter of dynamical system)
global N        % Number of segments
global FFTip    % True if simulation is for a follower force at the tip
global FFLength % True if simulation is for a follower force along the filament length
global Nf       % Number of filaments
global d        % Distance between filaments

% Simulation parameters (to edit)
% FIXED PARAMETERS
Nf = 2;
N = 20;
ndts = 1000;
% DO NOT ALTER

FFTip = true;
FFLength = false;

% ALTER THESE
f = fa;
d = da;

dt = 0.05;
Np = 20; % leave fixed
TotalSteps = 1000;

% END ALTER THESE
inputu = zeros(57,2);
pert = load("fixedperturbation.mat");
inputu(2:3:end, 2) = pert.pert;
inputu(2:3:end, 1) = inputu(2:3:end, 2);

out = initialvalueproblem2(f,reshape(inputu,3*(N-1),[]),dt,N,Nf,d,1000,FFTip,FFLength,0);

dataperiod = periodestimate(out(:, end-1), 0.05, false);

T = dataperiod(1);
start_dt = dataperiod(end-1);
T_dt = dataperiod(end) - start_dt;

Nphase = 16
u = zeros(TotalSteps, 3*(N-1)*Nf, Nphase);

parfor i=0:(Nphase-1)
    phase = i/Nphase
    
    end_dt = round(start_dt+phase*T_dt);
    
    phaseinput =  [out(start_dt, 1:3*(N-1)), out(end_dt, 3*(N-1)+1:end)]';

    u(:, :, i+1) = initialvalueproblem2(f,reshape(phaseinput,3*(N-1),[]),dt,N,Nf,d,TotalSteps,FFTip,FFLength,0);
end
save(sprintf('searchphase_f_%i_d_%i.mat',f,d), ...
    "f", "d", "u", "out");
end
