%-------------------------------------------------------------------------
%  advance x by ndts_ timesteps
%-------------------------------------------------------------------------
 function y = steporbit(ndts_,x)
   persistent dt
   global f
   global N
   global FFTip
   global FFLength
   global Nf
   global d

   if ndts_ ~= 1		% Set timestep size dt=T/ndts_
      dt = x(1) / ndts_ 	% If only doing one step to calc \dot{x},
   end				% then use previously set dt.

   a = x(2:end) ; %the initial Lie algebra set-up
   
   a = initialvalueproblem2(f,reshape(a,3*(N-1),[]),dt,N,Nf,d,ndts_,FFTip,FFLength,0);

   a = a(end,:)';
   y = zeros(size(x)) ;
   y(2:end) = a ;
 end 
 

