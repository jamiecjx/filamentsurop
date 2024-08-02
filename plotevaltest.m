[X,Y] = meshgrid(88:-1:66, 41:100);

surf(X, Y, real(evalarray(:, :, 1)))

set(gca, 'YDir','reverse')