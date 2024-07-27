function continuation_m(i)
    d = 89-i
    data = load(sprintf('f0_41_f1_50_d_%i_df_1_phase_5.000000e-01_ndts_400.mat', d))
    continuationarray = data.continuationarray;
    phase=data.phase;
    ndts=data.ndts;
    f0=50
    f1=100
    numf=50
    new_x = continuationarray(:, end);

    save("continuationinput.mat", "d", "f0", "f1", "numf", "new_x", "ndts", "phase") 

    continuation
end