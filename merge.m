% for d=66:88
%     w1 = load(sprintf('f0_41_f1_200_d_%i_df_1_phase_5.000000e-01_ndts_400.mat', d));
%     w2 = load(sprintf('f0_200_f1_300_d_%i_df_1_phase_5.000000e-01_ndts_400.mat', d));
%     delete(sprintf('f0_41_f1_200_d_%i_df_1_phase_5.000000e-01_ndts_400.mat', d));
%     delete(sprintf('f0_200_f1_300_d_%i_df_1_phase_5.000000e-01_ndts_400.mat', d));
%     cont1 = w1.continuationarray;
%     cont2 = w2.continuationarray;
%     df = 1;
%     f0 = 41;
%     f1 = 200;
%     ndts = 400;
%     phase = 0.5;
%     norm(cont1(:, end) - cont2(:, 1))
%     continuationarray = [cont1(:, 1:end-1), cont2];
%     save(sprintf('f0_%i_f1_%i_d_%i_df_%i_phase_%i_ndts_%i.mat',f0,f1,d,df,phase,ndts), ...
%     "continuationarray", "f0", "f1", "d", "df", "phase", "ndts");
% end

% for d=66:88
%     w1 = load(sprintf('f0_41_f1_200_d_%i_df_1_phase_5.000000e-01_ndts_400.mat', d));
%     delete(sprintf('f0_41_f1_200_d_%i_df_1_phase_5.000000e-01_ndts_400.mat', d));
%     cont1 = w1.continuationarray;
%     df = 1;
%     f0 = 41;
%     f1 = 300;
%     ndts = 400;
%     phase = 0.5;
%     continuationarray = w1.continuationarray;
%     save(sprintf('f0_%i_f1_%i_d_%i_df_%i_phase_%i_ndts_%i.mat',f0,f1,d,df,phase,ndts), ...
%     "continuationarray", "f0", "f1", "d", "df", "phase", "ndts");
% end



continuationarray = zeros(115, 260, 23);
for d=66:88
    w1 = load(sprintf('f0_41_f1_300_d_%i_df_1_phase_5.000000e-01_ndts_400.mat', d));
    continuationarray(:, :, 89-d) = w1.continuationarray;
end
d0=88
d1=66
f0=41
f1=300
df=1
dd=-1
phase=0.5
ndts=400
save(sprintf('2d_f0_%i_f1_%i_df_%i_d0_%i_d1_%i_dd_%i_phase_%i_ndts_%i.mat', ...
    f0, f1, df, d0, d1, dd, phase, ndts), "continuationarray", 'f0', 'f1', 'df', 'd0', 'd1', 'dd', 'phase', 'ndts')


% w1 = load('2d_f0_41_f1_100_df_1_d0_88_d1_66_dd_-1_phase_5.000000e-01_ndts_400.mat');
% w2 = load('2d_f0_100_f1_200_df_1_d0_88_d1_66_dd_-1_phase_5.000000e-01_ndts_400.mat');
% w3 = load('2d_f0_200_f1_300_df_1_d0_88_d1_66_dd_-1_phase_5.000000e-01_ndts_400.mat');
% continuationarray = zeros(115, 260, 23);
% continuationarray(:, 1:60, :) = w1;
% continuationarray(:, 61:160) = 