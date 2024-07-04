%*************************************************************************
% Openpipeflow.org.  If used in your work, please cite
% Willis, A. (2017) SoftwareX 6, 124-127.
% https://doi.org/10.1016/j.softx.2017.05.003 (open access)
%                                      Thanks in advance! Ashley 2019.
%*************************************************************************
% EXAMPLE: Periodic orbits of Lorenz system.
%- - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - -
%  Newton vector, current guess x_n:
%    x(1)     = T
%    x(2:end) = (x,y,z)
%
%  Function to minimise F(x):
%    F(1)     = 0.
%    F(2:end) = X(x,T) - x, difference between start and end point.
%
%  Newton update: x_{n+1} = x_n + dx.  Constraint on update dx:
%    dx . \dot{x} = 0, no update in direction of trajectory.
%
%*************************************************************************
%*************************************************************************
% PROGRAM MAIN
%*************************************************************************
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
% FIXED PARAMETERS
Nf = 2;
N = 20;
ndts = 200;

% DISTANCE APART

% d = 2*2.2*N;
% 
% 
% t0 = 9.65;
% x0 = jfnkstartantiphase;
% new_x = [t0, x0]'
% % FOLLOWER FORCE
% f = 100;
data = load("antiphasedatajfnk.mat")
d = data.d;
f = data.f;
new_x = data.new_x;

% DO NOT ALTER

FFTip = true;
FFLength = false;
mgmres  = 5;	% max GMRES iterations
nits    = 150 ;	    % max Newton iterations
rel_err = 1d-6;%1d-8;	% Relative error |F|/|x| %% i've changed this to be the error |F|: edited bottom of this script, and saveorbit

del     = -1d0 ;	% These rarely need changing for any problem
mndl    = 1d-20 ;
mxdl    = 1d+20 ;
gtol    = 1d-4;%1d-3 ;
epsJ    = 1d-6;%1d-5 ;

fixT = 0; % NOT FIXED POINTS

% PARAMETERS THAT ARE DETERMINED BY OLD PARAMETERS
n       = 3*(N-1)*Nf+1;	% Dimension of system, including unknown params


% plot initial guess
% plot_result

% scale parameters by |x| then call black box
ds = sqrt(dotprd(-1,new_x,new_x)) ; % had to change this: was overwriting d...
%% do we want rel_err?
tol  = rel_err * ds ;
del  = del     * ds ;
mndl = mndl    * ds ;
mxdl = mxdl    * ds ;

info = 0 ;
info = NewtonHook(@getrhs, @multJ, @multJp, @saveorbit, @dotprd, ...
               mgmres, n, gtol, tol, del, mndl, mxdl, nits, info) ;

% plot final solution
% plot_result
save("JFNKoutput.mat")

%*************************************************************************
% END PROGRAM MAIN
%*************************************************************************

function plot_result()

global new_x		% Current best x
global epsJ		% epsilon used in Jacobian approximation
global ndts		% Number of timesteps taken in period T
global fixT		% Fix T for equilibrium, rather than PO solution
global f		% Follower force (parameter of dynamical system)
global N
global Nf
global d
for i=Nf:-1:1
    v = reshape(new_x(2:end),3,[]);
    seed = [(i-1.5)*d;0;1.2];
    positions = GetPositionsFromEffectiveLieAlgebra(v(:,(N-1)*(i-1) + 1:(N-1)*i),seed);

    figure(2)
    PlotPositions(positions)
end
end


function x = GetPositionsFromEffectiveLieAlgebra(v,seed)

global N 

fil = Filament(N);

fil.X(:,1) = seed;

fac = 2^-0.5; % Equilibrium quaternion
Q_equilibrium = [fac*ones(N,1),zeros(N,1),-fac*ones(N,1),zeros(N,1)];

Qtemp = zeros(N,4);

Qtemp(1,:) = QuaternionProduct(qexp(zeros(3,1)),Q_equilibrium(1,:));

for i=2:N

    Qtemp(i,:) = QuaternionProduct(qexp(v(:,i-1)),Q_equilibrium(i,:));

end

fil.Q = [Qtemp, Qtemp];

fil.RobotArm;

x = fil.X;

end

function PlotPositions(x)

global N 

[s1,s2,s3] = sphere;

a = 1;

for j=1:N

    surf(x(1,j)+a*s1,x(2,j)+a*s2,x(3,j)+a*s3);
    
    hold on;

end

axis equal;

zlim([-1, 2.2*N+10]); xlim([-N-10,N+10]); ylim([-N-10,N+10]);

xlabel('X'); ylabel('Y'); zlabel('Z')

view([0 0 1]);
%view(0,0);

end