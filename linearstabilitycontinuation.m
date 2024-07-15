function linearstabilitycontinuation(file)
    global d
    global f
    global new_x
    
    data = load(file)
    
    % save("continuationinputantiphase.mat", "d", "f0", "f1", "numf", "new_x")
    continuationarray = data.continuationarray;
    
    d = data.d;
    f0 = data.f0;
    f1 = data.f1;
    df = data.df;
    numf = (f1-f0)/df;
    ndts = data.ndts;
    
    evalarray = zeros(numf+1, 6);
    evecarray = zeros(numf+1, 114, 6);
    
    parfor i=1:numf+1
        f = f0+df*(i-1);
        new_x = continuationarray(:, i);
        fprintf('floquet: starting iteration with f=%i\n',f) ;
        [eval, evec] = LinearStability(f, d, new_x(1), new_x(2:end), ndts);
        evalarray(i, :) = eval;
        evecarray(i, :, :) = evec;
    end
    
    save(sprintf('linearstability_f0_%i_f1_%i_d_%i_df_%i_phase_0.5.mat',f0,f1,d,df), ...
        "evalarray", "evecarray", "f0", "f1", "d", "df");
end