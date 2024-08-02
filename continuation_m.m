function continuation_m(i)
    d = 89-i
    data = load(sprintf('f0_41_f1_300_d_%i_df_1_phase_5.000000e-01_ndts_400.mat', d));
    continuationarray = data.continuationarray;
    phase=data.phase;
    ndts=data.ndts;
    f0=300;
    f1=400;
    numf=100;
    new_x = continuationarray(:, end);

    save(sprintf("continuationinput_%i.mat", i), "d", "f0", "f1", "numf", "new_x", "ndts", "phase") 

    continuation(sprintf("continuationinput_%i.mat", i))
    delete(sprintf("continuationinput_%i.mat", i))
end