function plot_solution(new_x)

		% Current best x
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

zlim([-1, 2.2*N+50]); xlim([-N-50,N+50]); ylim([-N-50,N+50]);

xlabel('X'); ylabel('Y'); zlabel('Z')

view([0 0 1]);
%view(0,0);

end