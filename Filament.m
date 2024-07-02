classdef Filament < handle
    
    % N.B. This is a handle class and does NOT have a copy constructor.
    % Assignment to a Filament type will cause both variables to be a
    % reference to the same memory location.
    
    properties
        
        Np; % Number of beads comprising the filament.
        
        Kb; % The elastic bending modulus of the filament.
        
        Kt; % The elastic twisting modulus of the filament.
        
        DeltaL; % The centreline distance between adjacent beads.
        
        StrainTwist; % The strain-twist vector encoding the intrinsic curvature and twist of the filament.
                        
        R; % The common radius of the constituent particles.
        
        X; % Positions of the particles comprising the filament.
        
        Q; % Orientation quaternions for the particles comprising the filament.
        
        U; % Lie algebra elements associated with the orientation of the constituent particles.
        
        V; % Velocities (both angular and translational) of the particles.
        
        F; % Forces and torques on the particles.
        
        Lambda; % The collection of Lagrange multipliers associated with the inextensibility of the filament.
        
        TetherLam; % The vector Lagrange multiplier associated with the tethering constraint.
        
        RotLam; % The vector Lagrange multiplier associated with the imposed rotation.
        
        Lmat; % The lower-triangular part of the Jacobian L-U decomposition.
        
        Umat; % The upper-triangular part of the Jacobian L-U decomposition.

        InitialForcesAndTorques;

        Xm1;

        Xm2;
        
    end
    
    methods
        
        function obj = Filament(varargin)
            
            if nargin==0
                
                % Empty constructor doesn't actually need to do anything,
                % it's only so that we can pre-allocate arrays to contain
                % Filaments.
                
            else
                
                N = varargin{1};
                
                obj.Np = N;
                
                obj.R = 1;

                obj.X = zeros(3,N);

                obj.Xm1 = zeros(3,N);
                
                obj.Xm2 = zeros(3,N);
                
                obj.Q = zeros(N,8);
                
                obj.U = zeros(9,N);
                
                obj.V = zeros(6,N);
                
                obj.F = zeros(6,N);
                
                obj.DeltaL = 2.2*obj.R;
                
                obj.Lambda = zeros(3,N-1);
                
                obj.TetherLam = zeros(3,1);
                
                obj.RotLam = zeros(3,1);

                obj.Kb = 18000;

                obj.Kt = obj.Kb;
                
            end
            
        end
        
        function RobotArm(obj)
            
            Xtemp = obj.X;
            
            Qtemp = obj.Q;
            
            dL = obj.DeltaL;
            
            for i=2:obj.Np
                
                Xtemp(:,i) = Xtemp(:,i-1) + ...
                    0.5*dL*(QuaternionRotation(Qtemp(i-1,1:4),[1;0;0]) + ...
                    QuaternionRotation(Qtemp(i,1:4),[1;0;0]));
                
            end
            
            obj.X = Xtemp;
            
        end
        
        function InitialSetup(obj,FirstBeadPos,StrainTwist,u)
            
            obj.X(:,1) = FirstBeadPos; % Because the filaments are tethered,
            % the position of the first bead is set here and never changed again.
            
            N = obj.Np;
  
            % All we need to restart from file is the quaternion described
            % by the 'effective Lie algebra' element

            u = reshape(u,3,[]); % = [u2,...,uN]

            fac = 2^-0.5; % Equilibrium quaternion
            Q_equilibrium = [fac*ones(N,1),zeros(N,1),-fac*ones(N,1),zeros(N,1)]; 

            Qtemp = zeros(N,4);

            Qtemp(1,:) = QuaternionProduct(qexp(zeros(3,1)),Q_equilibrium(1,:)); % u1 = [0;0;0] because of clamping

            for i=2:N

                Qtemp(i,:) = QuaternionProduct(qexp(u(:,i-1)),Q_equilibrium(i,:));

            end

            obj.Q = [Qtemp, Qtemp];

            obj.RobotArm;
            
            obj.StrainTwist = StrainTwist;

        end
        
        function InitialGuess(obj)
            
            Utemp = obj.U;
            
            Qtemp = obj.Q;
            
            for i=1:obj.Np
                
                % Guess Lie algebra elements, and hence quaternions, of all
                % particles.
                
                if i>1

                    Utemp(1:3,i) = 2*Utemp(4:6,i) - Utemp(7:9,i);
                    
                end
                
                Qtemp(i,1:4) = QuaternionProduct(qexp(Utemp(1:3,i)),Qtemp(i,5:8));
                
            end
            
            obj.U = Utemp;
            
            obj.Q = Qtemp;

            obj.Xm2 = obj.Xm1;

            obj.Xm1 = obj.X;
                        
        end
        
        function InternalForcesAndTorques(obj,f,FFTip,FFLength) % 
            
            N = obj.Np;
            
            ForcesAndTorques = zeros(6,N);
            
            Qtemp = obj.Q;
            
            Lam = obj.Lambda;
                        
            dL = obj.DeltaL;
            
            L = N * dL;

            ST = obj.StrainTwist;
            
            Bend = obj.Kb;
            
            Twist = obj.Kt;
            
            Tether = obj.TetherLam;
            
            Rot = obj.RotLam;
            
            for i=1:obj.Np-1
                
                % Constraint forces and torques
                
                ForcesAndTorques(1:3,i) = ForcesAndTorques(1:3,i) - Lam(:,i);
                
                ForcesAndTorques(1:3,i+1) = ForcesAndTorques(1:3,i+1) + Lam(:,i);
                
                t = QuaternionRotation(Qtemp(i,1:4),[1;0;0]);
                
                ForcesAndTorques(4:6,i) = ForcesAndTorques(4:6,i) - ...
                    0.5*dL*[t(2)*Lam(3,i) - t(3)*Lam(2,i);t(3)*Lam(1,i) - t(1)*Lam(3,i);t(1)*Lam(2,i) - t(2)*Lam(1,i)];
                
                t = QuaternionRotation(Qtemp(i+1,1:4),[1;0;0]);
                
                ForcesAndTorques(4:6,i+1) = ForcesAndTorques(4:6,i+1) - ...
                    0.5*dL*[t(2)*Lam(3,i) - t(3)*Lam(2,i);t(3)*Lam(1,i) - t(1)*Lam(3,i);t(1)*Lam(2,i) - t(2)*Lam(1,i)];
                
                % Elastic torques
                
                q = MidpointQ(Qtemp(i,1:4),Qtemp(i+1,1:4));
                
                dqds = (Qtemp(i+1,1:4) - Qtemp(i,1:4))/dL;
                
                b = 2*QuaternionProduct([q(1),-q(2),-q(3),-q(4)],dqds);
                b = b(2:4);
                
                M = QuaternionRotation(q,[Twist * (b(1) - ST(3)); ...
                    Bend * (b(2) - ST(1)); Bend * (b(3) - ST(2))]);
                
                ForcesAndTorques(4:6,i) = ForcesAndTorques(4:6,i) + M;
                
                ForcesAndTorques(4:6,i+1) = ForcesAndTorques(4:6,i+1) - M;
                
            end

            % Follower-force
            FollowerForce = -f*Bend/(L*L);

            if FFTip
                
                tN = QuaternionRotation(Qtemp(N,1:4),[1;0;0]);

                ForcesAndTorques(1:3,N) = ForcesAndTorques(1:3,N) + FollowerForce*tN;
            
            elseif FFLength

                for i=1:N

                    t = QuaternionRotation(Qtemp(i,1:4),[1;0;0]);

                    const = 1;%(i-1)/(N-1);

                    ForcesAndTorques(1:3,i) = ForcesAndTorques(1:3,i) + const*FollowerForce*t;

                end
            
            end

            
            % Tethering force and torque
            
            u = obj.U(1:3,1);
            
            ForcesAndTorques(1:3,1) = ForcesAndTorques(1:3,1) + Tether;
            
            theta = (u(1)*u(1) + u(2)*u(2) + u(3)*u(3))^0.5;
            
            if theta<10^-6
                fac = -1/12;
            else
                fac = (0.5*theta*cot(0.5*theta) - 1)/(theta^2);
            end
            
            TetherTorque = Rot - 0.5*cross(Rot,u) - fac*cross(cross(Rot,u),u);
            
            ForcesAndTorques(4:6,1) = ForcesAndTorques(4:6,1) + TetherTorque;

            obj.F = ForcesAndTorques;
            
        end

        function Y = GetPositions(obj)

            Y = obj.X;

        end
        
        function com = CentreOfMass(obj)
            
            com = mean(obj.X,2);
            
        end
        
        function DecomposeJacobian(obj,dt,mu,FFTip,FFLength,f)
            
            [obj.Lmat,obj.Umat] = lu(ConstructApproxJ(obj,dt,mu,FFTip,FFLength,f));
            
        end
        
        function out = InvertLocalBlock(obj,v)
            
            out = obj.Lmat\v;
            
            out = obj.Umat\out;
            
        end
        
        function ApplyUpdate(obj,u)
            
            N = obj.Np;
            
            obj.TetherLam = obj.TetherLam + u(1:3);
            
            obj.RotLam = obj.RotLam + u(4:6);
            
            Utemp = obj.U;
            
            Qtemp = obj.Q;
            
            Lam = obj.Lambda;
            
            for i=2:N
                
                Utemp(1:3,i) = Utemp(1:3,i) + u(3*i+1:3*(i+1));
                
                Qtemp(i,1:4) = QuaternionProduct(qexp(Utemp(1:3,i)),Qtemp(i,5:8));
                
            end
            
            for i=1:N-1
                
                Lam(:,i) = Lam(:,i) + u(3*(N+i)+1:3*(N+i+1));
                
            end
            
            obj.U = Utemp;
            
            obj.Q = Qtemp;
            
            obj.Lambda = Lam;
            
        end
        
        function Lam = EndOfStepUpdate(obj)
            
            Utemp = obj.U;
            
            Qtemp = obj.Q;
            
            LamTemp = obj.Lambda;

            Lam = LamTemp(:,end);
            
            Utemp(7:9,:) = Utemp(4:6,:);
            
            Utemp(4:6,:) = Utemp(1:3,:);
            
            for i=1:obj.Np
                
                Qtemp(i,5:8) = QuaternionProduct(qexp(Utemp(1:3,i)),Qtemp(i,5:8));
                
            end
            
            obj.U = Utemp;
            
            obj.Q = Qtemp;

            % Vtemp = obj.V;
            % Xtemp = obj.X;
            % temp = zeros(1,obj.Np);

            % for i=1:obj.Np
            % 
            %     ti = QuaternionRotation(Qtemp(i,1:4),[1;0;0]);
            % 
            %     % omega = norm(cross(Vtemp(4:6,i),ti))/norm(cross([0;0;1],ti));
            %     if abs(ti(1))>abs(ti(2))
            % 
            %         omega = Vtemp(6,i) - Vtemp(4,i)*ti(3)/ti(1);
            % 
            %     else
            % 
            %         omega = Vtemp(6,i) - Vtemp(5,i)*ti(3)/ti(2);
            % 
            %     end
            % 
            %     temp(i) = norm(Vtemp(1:3,i) - cross(omega*[0;0;1],Xtemp(:,i)));
            % 
            % end
            % temp
        end
            
        function x = GetEffectiveLieAlgebraUpdate(obj)
        
            Qn = obj.Q(:,1:4);
            N = obj.Np;
            x = zeros(3,N-1);
            fac = 2^-0.5;
            Q_equilibrium_adjoint_i = [fac,0,fac,0];
          
            % Want to find the 'effective' Lie algebra element that gives our
            % segment positions w.r.t the vertical equilibrium
            for i=2:N % ignoring first entry ([0,0,0] because of clamping condition)
            
                dot_product = QuaternionProduct(Qn(i,:),Q_equilibrium_adjoint_i);
                
                sinmodv2 = norm(dot_product(2:4));
            
                if sinmodv2 ~= 0
                        
                    xi = dot_product(2:4)/sinmodv2;
            
                else
            
                    xi = zeros(1,3);
                
                end
            
                if dot_product(1) > 1

                    modvtemp = 0;
                    
                else

                    modvtemp = 2*acos(dot_product(1));
                
                end

                x(:,i-1) = xi*modvtemp;
            
            end

            x = reshape(x,3*(N-1),[]);

        end
       
        function [Vtemp,Lam] = GetVelocity(obj)

            Vtemp = obj.V;
            Lam = obj.Lambda;

        end
    end
    
end

