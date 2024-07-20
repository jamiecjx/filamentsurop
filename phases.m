function output = phases(u)
    output = zeros(1, size(u, 3));
    for i=1:size(u,3)
        [y1, y2] = phaseestimate(u(:,:,i))
        j = size(y2, 2);
        while y2(j) >= y1(end)
            j = j-1;
        end
        output(i) = (y1(end-1) - y2(j)) / (y1(end-1) - y1(end));
    end
end