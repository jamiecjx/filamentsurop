function [evalue,evector] = LinearStability(f,d,T,xstar,ndts)

% inputs: f = dimensionless follower force
%         d = distance between filaments
%         T = period
%         xstar = equilibrium configuration
%         ndts = total timesteps per iteration

% Define Simulation type
FFTip=true;
FFLength=false;

% Set-up linear stability parameters
epsilon = 1e-3; % size of perturbation

dt = T/ndts;

tol = 1e-8;%1e-6; % error tolerance
difference = 1;
Number_Evals_We_Want = 6;

% Set-up Filament code parameters
N = 20; Nf=2;

V = zeros(3*(N-1)*Nf,1);
H = zeros(3*(N-1)*Nf,3*(N-1));
evalue = zeros(1,Number_Evals_We_Want);
oldevals = 100*ones(Number_Evals_We_Want,1);

% Setting up initial U perturbation
b = zeros(3,Nf*(N-1));
b(2,:) = randn(1,Nf*(N-1));
% b(1:2,:) = randn(2,N-1);
% b = randn(57,1)
b = reshape(b,3*(N-1)*Nf,[]);
q1 = b/norm(b);

V(:,1) = q1;

iter = 0;
while difference>tol

    iter = iter+1;
    
    v = initialvalueproblem2(f,reshape(xstar + epsilon*V(:,iter),3*(N-1),[]),dt,N,Nf,d,ndts,FFTip,FFLength,0);

    v = v(end,:)';

    v = reshape(v,Nf*3*(N-1),[]);
    
    v = (v-xstar)/epsilon;

    for i=1:iter

        H(i,iter) = V(:,i)'*v;

        v = v - H(i,iter) * V(:,i);

    end

    H(iter+1,iter) = norm(v);
    
    if iter >= Number_Evals_We_Want
        
        evalues = eig(H(1:iter,1:iter));

        fprintf('After log and sorting:')

        evalues(1:iter) = log(evalues(1:iter))/T

        evalues = sort(evalues,'descend','ComparisonMethod','real');

        difference = norm(evalues(1:Number_Evals_We_Want)-oldevals)
    
        oldevals = evalues(1:Number_Evals_We_Want);

        if H(iter+1,iter) < 1e-14

            disp('Arnoldi converged! Eigenvalues of H = eigenvalues of Jacobian')
            
            evalue(1:Number_Evals_We_Want) = evalues(1:Number_Evals_We_Want);

            [evec,eval] = eig(H(1:iter,1:iter));

            [~,ideval] = sort(diag(eval),'descend');

            evec_sorted = evec(:,ideval);

            evector = V*evec_sorted;

            evector = evector(:,1:Number_Evals_We_Want);
 
            return
        end

    end

    V(:,iter+1) = v/H(iter+1,iter);

end
  
k = min(Number_Evals_We_Want,length(evalues));
evalue(1:k) = evalues(1:k);

[evec,eval] = eig(H(1:iter,1:iter));
[~,ideval] = sort(diag(eval),'descend');

evec_sorted = evec(:,ideval);

evector = V(:,1:iter)*evec_sorted;

evector = evector(:,1:Number_Evals_We_Want);

end

% ---------------------------------------------------------------------
% Local functions
% ---------------------------------------------------------------------

