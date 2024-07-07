%-------------------------------------------------------------------------
%  estimate period of solution based off angle of topmost point of filament
%-------------------------------------------------------------------------
 function p = periodestimate(x, dt, plotgraph)
     if plotgraph
         plot(x)
     end
     y = [];
     for i = 2:length(x)-1
         if x(i-1) < x(i) && x(i) > x(i+1)
            y = [y, i];
         end
     end
     differences = diff(y);
     y
     p = [differences(end) * dt, y];
 end

