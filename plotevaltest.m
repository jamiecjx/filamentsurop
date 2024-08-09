[X,Y] = meshgrid(88:-1:66, 301:400);
hold on
surf(X, Y, real(evalarray(:, :, 1)))
surf(X, Y, real(evalarray(:, :, 2)))

set(gca, 'YDir','reverse')