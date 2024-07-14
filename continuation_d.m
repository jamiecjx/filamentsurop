
function continuation_d(i)
    global d
    global f
    global new_x
    global ndts
    data = load("continuation_d_input.mat")
    
    % save("continuation_d_input.mat", "d0", "d1", "f0", "f1", "numf", "numd", "continuationarray", "ndts", "phase")
    
    d0 = data.d0;
    d1 = data.d1;
    f0 = data.f0;
    f1 = data.f1;
    numf = data.numf;
    numd = data.numd;
    continuationarray = data.continuationarray;
    ndts = data.ndts;
    phase = data.phase;
    
    df = (f1-f0)/numf;
    dd = (d1-d0)/numd;
    
    continuationarray3d = zeros(size(new_x, 1), numf+1, numd+1);
    M = numf+1
    L = numd+1
    
    f = f0+df*(i-1);
    new_x = continuationarray(:,i);
    for j=1:L
        d = d0+dd*(j-1);
        fprintf('continuation: starting iteration with f=%i, d=%i',f,d) ;
        JFNK
        continuationarray3d(:,i,j) = new_x;
    end
    
    
    save(sprintf('f0_%i_f1_%i_d0_%i_d1_%i_df_%i_dd_%i_phase_%i_ndts_%i.mat',f0,f1,d0,d1,df,dd,phase,ndts), ...
        "continuationarray3d", "f0", "f1", "d", "df", "phase", "ndts");
end