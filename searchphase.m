function u = searchphase(fa, da, Nphase, TotalSteps, dt)

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
ndts = 200;
% DO NOT ALTER

FFTip = true;
FFLength = false;

% ALTER THESE
f = fa;
d = da;


Np = 20; % leave fixed

% END ALTER THESE
inputu = zeros(57,2);
pert = load("fixedperturbation.mat");
inputu(2:3:end, 2) = pert.pert;
inputu(2:3:end, 1) = inputu(2:3:end, 2);

out = initialvalueproblem2(f,reshape(inputu,3*(N-1),[]),dt,N,Nf,d,TotalSteps,FFTip, FFLength,0);

dataperiod = periodestimate(out(:, end-1), dt, false);

T = dataperiod(1);
start_dt = dataperiod(end-1);
T_dt = dataperiod(end) - start_dt;


un1 = zeros(TotalSteps, Nphase);
un2 = zeros(TotalSteps, Nphase);
ut = zeros(3*(N-1)*Nf, Nphase);
parfor i=0:(Nphase-1)
    phase = i/Nphase
    
    end_dt = round(start_dt+phase*T_dt);
    
    phaseinput =  [out(start_dt, 1:3*(N-1)), out(end_dt, 3*(N-1)+1:end)]';

    ivp = initialvalueproblem2(f,reshape(phaseinput,3*(N-1),[]),dt,N,Nf,d,TotalSteps,FFTip,FFLength,0);
    un1(:, i+ 1) = ivp(:, 113);
    un2(:, i+ 1) = ivp(:, 56);
    ut(:, i+ 1) = ivp(end, :);
end
save(sprintf('searchphase_f_%i_d_%i_Nphase_%i_ts_%i_ndts_%i_dt_%i.mat',f,d,Nphase,TotalSteps,ndts,dt), ...
    "f", "d", "un1", "un2", "ut", "out");
end
