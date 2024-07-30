%-------------------------------------------------------------------------
%  estimate phase of solution based off angle of topmost point of filament
%-------------------------------------------------------------------------
function  phaseestimate(un1, un2)
     function y=phasepoints(x)
         y = [];
         length(x);
         for i = 2:length(x)-1
             if x(i-1) < x(i) && x(i) > x(i+1)
                y = [y, i];
             end
         end
     end
        % hardcoded because lazy
        y1 = phasepoints(un1);
        y2 = phasepoints(un2);
        % hold on
        % plot(u(:, 113))
        % plot(u(:, 56))
        % hold off
        y = diff(sort([y1, y2]))
        L = size(y, 2)
        scatter(1:L, y)
 end

% latest result is converging to in phase.