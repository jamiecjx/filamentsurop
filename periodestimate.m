%-------------------------------------------------------------------------
%  estimate period of solution based off angle of topmost point of filament
%-------------------------------------------------------------------------
 function p = periodestimate(x, dt, plotgraph)
     if plotgraph
         plot(x)
     end
     y = [];
     length(x);
     for i = 2:length(x)-1
         if x(i-1) < x(i) && x(i) > x(i+1)
            y = [y, i];
         end
     end
     differences = diff(y);
     k = 80
     p = [sum(differences(end-k+1:end))/k * dt, y];
 end

% hold on
% periodestimate(out(:, 113), 0.0039, true)
% periodestimate(out(:, 56), 0.0039, true)
% hold off