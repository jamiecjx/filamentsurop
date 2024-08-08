function linearstabilitycontinuationd(i, file)
    data = load(file)

    continuationarray = data.continuationarray;
  
    d0 = data.d0;
    d1 = data.d1;
    f0 = data.f0;
    f1 = data.f1;
    df = data.df;
    dd = data.dd;
    phase = data.phase;

    numd = (d1-d0)/dd;
    L = numd+1;
    numf = (f1-f0)/df;
    ndts = data.ndts;
    
    evalarray = zeros(numd+1, 6);

    f = f0+df*(i-1);

    for j=1:L
        d_ = d0+dd*(j-1);
        new_xx = continuationarray(:, i, j);
        fprintf('floquet: starting iteration with f=%i, d=%i\n',f, d_) ;
        [eval, ~] = LinearStability(f, d_, new_xx(1), new_xx(2:end), ndts);
        evalarray(j, :) = eval;
    end

    save(sprintf('ls_f_%i_df_%i_d0_%i_d1_%i_dd_%i_phase_%i_ndts_%i.mat',f,df,d0,d1,dd,phase,ndts), ...
        "evalarray", "f0", "f1", "d0", "d1", "dd", "ndts", "phase", "df");
end
