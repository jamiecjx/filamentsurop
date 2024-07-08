function new_x = createphase(fa, da, phase)

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
phase
% END ALTER THESE

u = zeros(57,2);
pert = load("fixedperturbation.mat");
u(2:3:end, 2) = pert.pert;
u(2:3:end, 1) = u(2:3:end, 2);

out = initialvalueproblem2(f,reshape(u,3*(N-1),[]),4/f,N,Nf,d,500,FFTip,FFLength,0);

dataperiod = periodestimate(out(:, end-1), 0.05, true);

T = dataperiod(1);
start_dt = dataperiod(end-1);
T_dt = dataperiod(end) - start_dt;

end_dt = round(start_dt+phase*T_dt)

new_x =  [T, out(start_dt, 1:3*(N-1)), out(end_dt, 3*(N-1)+1:end)]';

JFNK

save(sprintf('jfnk_f_%i_d_%i_phase_%i.mat',f,d,phase), ...
    "f", "d", "new_x", "phase");
end
