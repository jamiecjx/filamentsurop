function [C] = J_block_C(Fibre,fac,FFTip,FFLength,dt,mu,a,f)

Np = Fibre.Np;

Q = Fibre.Q;

C = zeros(3*(Np-1),3*Np);

C(1:3,1:6) = fac*[rcross(QuaternionRotation(Q(1,1:4),[1;0;0])),rcross(QuaternionRotation(Q(2,1:4),[1;0;0]))];

for k=3:Np
    
    i = 3*(k-2);
    
    C(i+1:i+3,1:i+3) = C(i-2:i,1:i+3);
    
    C(i+1:i+3,i+1:i+6) = C(i+1:i+3,i+1:i+6) + fac*[rcross(QuaternionRotation(Q(k-1,1:4),[1;0;0])),rcross(QuaternionRotation(Q(k,1:4),[1;0;0]))];
    
end

if FFTip

elseif FFLength
    
    FollowerForce = -f*Fibre.Kb/(Fibre.DeltaL*Np)^2;

    const = dt*FollowerForce/(9*pi*mu*a);

    for k=2:Np
        tk = QuaternionRotation(Q(k,1:4),[1;0;0]);
        C(3*k-5:3*(k-1),3*k-2:3*k) = C(3*k-5:3*(k-1),3*k-2:3*k) + const*rcross(tk);
    end

end

end

