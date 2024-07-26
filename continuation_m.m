function continuation_m(i)
    data = load('f_41_d0_88_d1_66_df_1_dd_-1_phase_5.000000e-01_ndts_400.mat')
    continuationarray3d = data.continuationarray3d
    phase=data.phase
    ndts=data.ndts
    f0=41
    f1=50
    d=89-i
    numf=9
    new_x = continuationarray3d(:, i);

    save("continuationinput.mat", "d", "f0", "f1", "numf", "new_x", "ndts", "phase") 

    continuation
end