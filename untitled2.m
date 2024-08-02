for d=66:88
    data = load(sprintf('f0_50_f1_100_d_%i_df_1_phase_5.000000e-01_ndts_400.mat', d));
    continuationarray = data.continuationarray
    t = 51;
    while norm(continuationarray(:, t)) == 0
        t = t - 1;
    end
    t
end
