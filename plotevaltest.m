[X,Y] = meshgrid(88:-1:66, 41:400);
hold on
surf(X, Y, real(evalarray(:, :, 1)))
surf(X, Y, real(evalarray(:, :, 2)))

set(gca, 'YDir','reverse')


% data = load('ls_f0_41_f1_300_df_1_d0_88_d1_66_dd_-1_phase_5.000000e-01_ndts_400.mat')
% we1 = data.evalarray;
% data2 = load('ls_f0_301_f1_400_df_1_d0_88_d1_66_dd_-1_phase_5.000000e-01_ndts_400.mat')
% 
% we2 = data2.evalarray;
% size(we1)
% size(we2)
% evalarray = cat(1, we1, we2);
% 
% d0 = 88
% d1 = 66
% dd = -1
% df = 1
% f0 = 41
% f1=400
% ndts=400
% phase=0.5
% 
% save(sprintf('ls_f0_%i_f1_%i_df_%i_d0_%i_d1_%i_dd_%i_phase_%i_ndts_%i.mat',f0,f1,df,d0,d1,dd,phase,ndts), ...
%         "evalarray", "f0", "f1", "d0", "d1", "dd", "ndts", "phase", "df");