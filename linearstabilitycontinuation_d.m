function linearstabilitycontinuation_d(file)
    global d
    global f
    global new_x
    
    data = load(file)
    
    % save("continuationinputantiphase.mat", "d", "f0", "f1", "numf", "new_x")
    continuationarray3d = data.continuationarray3d;
    
    d0 = data.d0;
    d1 = data.d1;
    f = data.f;
    dd = data.dd;
    numd = (d1-d0)/dd;
    ndts = data.ndts;
    
    evalarray = zeros(numd+1, 6);
    evecarray = zeros(numd+1, 114, 6);
    
    parfor i=1:numd+1
        d = d0+dd*(i-1);
        new_x = continuationarray3d(:, i);
        fprintf('floquet: starting iteration with d=%i\n',d) ;
        [eval, evec] = LinearStability(f, d, new_x(1), new_x(2:end), ndts);
        evalarray(i, :) = eval;
        evecarray(i, :, :) = evec;
    end
    
    save(sprintf('linearstability_d_d0_%i_d1_%i_f_%i_dd_%i_phase_0.5.mat',d0,d1,f,dd), ...
        "evalarray", "evecarray", "d0", "d1", "f", "dd");
end