function continuation_m(i)
    d = 89-i
    data = load(sprintf('f0_200_f1_300_d_%i_df_1_phase_0_ndts_400.mat', d));
    continuationarray = data.continuationarray;
    phase=data.phase;
    ndts=data.ndts;
    f0=250;
    f1=300;
    numf=50;
    new_x = continuationarray(:, 51);

    save(sprintf("continuationinput_%i.mat", i), "d", "f0", "f1", "numf", "new_x", "ndts", "phase") 

    continuation(sprintf("continuationinput_%i.mat", i))
    delete(sprintf("continuationinput_%i.mat", i))
end

% function continuation_m(i)
%     d = 89-i
%     data = load('f_41_d0_88_d1_66_df_0_dd_-1_phase_0_ndts_400.mat');
%     continuationarray3d = data.continuationarray3d;
%     phase=data.phase;
%     ndts=data.ndts;
%     f0=41;
%     f1=50;
%     numf=9;
%     new_x = continuationarray3d(:, i);
% 
% 
%     save(sprintf("continuationinput_%i.mat", i), "d", "f0", "f1", "numf", "new_x", "ndts", "phase") 
% 
%     continuation(sprintf("continuationinput_%i.mat", i))
%     delete(sprintf("continuationinput_%i.mat", i))
% end