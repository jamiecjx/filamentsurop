% for i=40:50
%     f=i
%     d0=88
%     d1=66
%     dd=-1
%     phase=0.5
%     ndts=200
%     data = load(sprintf('f_%i_d0_%i_d1_%i_df_1_dd_-1_phase_%i_ndts_%i.mat',f,d0,d1,phase,ndts))
%     continuationarray=data.continuationarray3d
%     save(sprintf('f_%i_d0_%i_d1_%i_dd_%i_phase_%i_ndts_%i.mat',f,d0,d1,dd,phase,ndts), ...
%         "continuationarray", "f", "d0", "d1", "dd", "phase", "ndts");
% end

% data=load("linearstability_f0_200_f1_450_d_88_df_1_phase_0.5.mat")
% weval = data.evalarray;
% wevec = data.evecarray;
% 
% data2=load("linearstability_f0_450_f1_500_d_88_df_1_phase_0.5.mat")
% 
% evalarray = data2.evalarray;
% evecarray = data2.evecarray;
% 
% evalarray = cat(1, weval, evalarray(2:end, :));
% evecarray = cat(1, wevec, evecarray(2:end, :, :));
% 
% f0=200;
% f1=500;
% df=1;
% d=88;
% 
% save(sprintf('linearstability_f0_%i_f1_%i_d_%i_df_%i_phase_0.5.mat',f0,f1,d,df), ...
%         "evalarray", "evecarray", "f0", "f1", "d", "df");


% save(sprintf('linearstability_f0_%i_f1_%i_d_%i_df_%i_phase_0.5.mat',f0,f1,d,df), ...
%         "evalarray", "evecarray", "f0", "f1", "d", "df");

% continuationarray2d = zeros(115, 23, 10);
% for i=41:50
%     f=i
%     d0=88
%     d1=66
%     dd=-1
%     phase=0.5
%     ndts=200
%     data = load(sprintf('f_%i_d0_%i_d1_%i_dd_-1_phase_%i_ndts_%i.mat',f,d0,d1,phase,ndts))
%     continuationarray2d(:, :, i-40) = data.continuationarray
% end
% 
% f0=41
% f1=50
% df=1
% d0=88
% d1=66
% dd=-1
% phase=0.5
% ndts=200
% save(sprintf('2d_f0_%i_f1_%i_df_%i_d0_%i_d1_%i_dd_%i_phase_%i_ndts_%i.mat',f0, f1, df,d0,d1,dd,phase,ndts), ...
%          "continuationarray2d", "f0", "f1", "df", "d0", "d1", "dd", "phase", "ndts");



% hold on
% for f=41:50
%     data=load(sprintf('linearstability_d_d0_88_d1_66_f_%i_dd_-1_phase_5.000000e-01_ndts_200.mat', f));
%     evalarray= data.evalarray;
%     plot(real(evalarray(:, 1:2)), 'DisplayName',string(f),'Color',[(f-41)/9, 1, 0])
% end



data = load('f0_200_f1_500_d_88_df_1_phase_5.000000e-01_ndts_400.mat');
wh = data.continuationarray;
data = load('f0_500_f1_550_d_88_df_1_phase_5.000000e-01_ndts_400.mat');
continuationarray = data.continuationarray;
disp(norm(wh(:, end) - continuationarray(:, 1)))
continuationarray = [wh, continuationarray(:, 2:end)];

f0=200
f1=550
d=88
df=1
phase=0.5
ndts=400
save(sprintf('f0_%i_f1_%i_d_%i_df_%i_phase_%i_ndts_%i.mat',f0,f1,d,df,phase,ndts), ...
    "continuationarray", "f0", "f1", "d", "df", "phase", "ndts");