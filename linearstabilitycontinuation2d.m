function linearstabilitycontinuation2d(file)
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
    
    evalarray = zeros(numf+1, numd+1, 6);
    parpool
    parfor i=1:numf+1
        for j=1:L
            f_ = f0+df*(i-1);
            d_ = d0+dd*(j-1);
            new_xx = continuationarray(:, i, j);
            fprintf('floquet: starting iteration with f=%i, d=%i\n',f_, d_) ;
            [eval, ~] = LinearStability(f_, d_, new_xx(1), new_xx(2:end), ndts);
            evalarray(i, j, :) = eval;
        end
    end
    
    save(sprintf('ls_f0_%i_f1_%i_df_%i_d0_%i_d1_%i_dd_%i_phase_%i_ndts_%i.mat',f0,f1,df,d0,d1,dd,phase,ndts), ...
        "evalarray", "f0", "f1", "d0", "d1", "dd", "ndts", "phase", "df");
end