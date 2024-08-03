% function [EffLieAlgebra] = main2(f,u,dt,Np,Nf,d,TotalSteps,FFTip,FFLength,vid)
function [EffLieAlgebra] = main2

% Progressing simulation forward by 1 timestep

% inputs: f = dimensionless follower force
%         u = Lie algebra element defining quaternion w.r.t vertical
%         equilibrium (u2,...,uN)
%         dt = timestep
%         Np = number of segments
%         Nf = number of filaments
%         d = distance between filaments
%         TotalSteps = number of timesteps
%         FFTip = 1 if follower force at the tip, 0 otherwise
%         FFLength = 1 if follower force along the length, 0 otherwise
%         vid = 1 to plot filament and save video, 0 otherwise
global a
global lol
% Variables to get started:

% 
% u = reshape(a(2:end), 57, [])
%u = reshape(lol, 57, []);
% data = load('store.mat')
% u = data.storeyooo;
% u = a;
% u = u(:, end);
% u = u(2:end);
% u = reshape(u, 57, [])

dt = 0.05;
Np = 40; % leave fixed
TotalSteps = 100000;
FFTip = 1; % leave fixed
FFLength = 0; % leave fixed
vid = 1;
Nf = 2;%2; % Number of filaments.
mu = 1; % Fluid viscosity. leave fixed
L = 2.2*Np;
f = 46;
d=44;

data = load("fixedperturbation40.mat")
u = zeros(3*Np-3, 2);
u(2:3:end, 1) = data.pert;
u(2:3:end, 2) = -u(2:3:end, 1)+ 0.01*rand()




if vid
    VideoName = sprintf('f_equals_%i_Np_%i_d_%i_video.avi',f,Np,d);
    try
        video = VideoWriter(VideoName, 'MPEG-4');
    catch
        video = VideoWriter(VideoName);
    end
    open(video);
end

for i=Nf:-1:1
    Filaments(i) = Filament(Np);
    Filaments(i).InitialSetup([(i-1.5)*d;0;1.2],[0,0,0],u(:,i));
    N(i) = 6*Filaments(i).Np;
end

N = [0,cumsum(N)];
Nbroy = N(end);

for Steps=1:TotalSteps

    Cmat = zeros(Nbroy,20); Dmat = zeros(Nbroy,20);
    
    for i=1:Nf
        
        Filaments(i).InitialGuess;
        Filaments(i).RobotArm;
        Filaments(i).InternalForcesAndTorques(f,FFTip,FFLength);
        
    end

    StericForces(Filaments);
    
    RPY(Filaments,mu)
    
    [Check,Error] = ConstraintCheck(Filaments,dt,Steps);
    
    for i=1:Nf
        
        Filaments(i).DecomposeJacobian(dt,mu,FFTip,FFLength,f);
        
    end
    
    BroydenIter = 0;
    
    while Check==1
        
        Update = ApplyInverseJacobian(-Error,Filaments,Cmat,Dmat,BroydenIter);
    
        for i=1:Nf
            
            Filaments(i).ApplyUpdate(Update(N(i)+1:N(i+1)));
            Filaments(i).RobotArm;
            Filaments(i).InternalForcesAndTorques(f,FFTip,FFLength);
            
        end
        
        StericForces(Filaments);
        
        RPY(Filaments,mu)
        
        [Check,NewError] = ConstraintCheck(Filaments,dt,Steps);
        
        Svec = ApplyInverseJacobian(-NewError,Filaments,Cmat,Dmat,BroydenIter);
        
        Diff = NewError - Error;
        
        BroydenIter = BroydenIter + 1;
        
        fac = norm(Diff);
        
        Cmat(:,BroydenIter) = Svec/fac;
        
        Dmat(:,BroydenIter) = Diff/fac;
        
        Error = NewError;
     
    end

    for i=1:Nf
        
        Lambda = Filaments(i).EndOfStepUpdate;
        Y(3*Steps-2:3*Steps,:) = Filaments(i).GetPositions;
        EffLieAlgebra(Steps,3*(Np-1)*(i-1)+1:3*(Np-1)*i) = Filaments(i).GetEffectiveLieAlgebraUpdate;

    end

    if mod(Steps, 50)==1
        fprintf('Step %i/%i required %i Broyden iterations.\n',Steps,TotalSteps,BroydenIter);
        save(sprintf('ivp_f_%i_d_%i',f,d), "EffLieAlgebra", "f", "d");
    end

    if vid && mod(Steps,10)==1
        MakePlot(Nf,Np,Filaments,video)
    end

end

if vid
    close(video)
end
save(sprintf('ivp_f_%i_d_%i_Np_%i_st_%i',f,d,Np, TotalSteps), "EffLieAlgebra", "f", "d");
end

%%%%%% LOCAL FUNCTIONS %%%%%%

function MakePlot(Nf,Np,Filaments,video)

% Plotting
    figure(2)
[s1,s2,s3] = sphere;

for i=1:Nf
    
    MyColour = ((i-1)/Nf)*[1,1,1];
    
    a = Filaments(i).R;
    
    x = Filaments(i).X;
    
    for j=1:Np
        
        surf(x(1,j)+a*s1,x(2,j)+a*s2,x(3,j)+a*s3,'FaceColor',MyColour);
        hold on;
        
    end

end
 
axis equal;

L = 2.2*Np;

zlim([-1, L+10]); xlim([-70,70]); ylim([-70,70]); 

xlabel('X'); ylabel('Y'); zlabel('Z');

view(0,0);
% view([0 0 1]);

CurrFrame = getframe(gcf);
writeVideo(video,CurrFrame);

pause(0.1)
hold off;
     
end

