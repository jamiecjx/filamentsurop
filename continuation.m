global d
global f
global new_x
global ndts
data = load("continuationinput.mat")

% save("continuationinput.mat", "d", "f0", "f1", "numf", "new_x", "ndts", "phase")



% d = 88;
% f0 = 100;
% f1 = 300;
% numf = 200;
% new_x = continuationarray(:,1);

d = data.d;
f0 = data.f0;
f1 = data.f1;
numf = data.numf;
new_x = data.new_x;
ndts = data.ndts;
phase = data.phase;

df = (f1-f0)/numf;

continuationarray = zeros(size(new_x, 1), numf+1);
continuationarray(:,1) = new_x;
sprintf('f0_%i_f1_%i_d_%i_df_%i_phase_0.mat',f0,f1,d,df)
for i=1:numf+1
    save(sprintf('progress_f0_%i_f1_%i_d_%i_df_%i_phase_%i_ndts_%i.mat',f0,f1,d,df,phase,ndts), ...
    "i");

    f = f0+df*(i-1);
    fprintf('continuation: starting iteration with f=%i',f) ;
    JFNK
    continuationarray(:,i) = new_x;
    save(sprintf('f0_%i_f1_%i_d_%i_df_%i_phase_%i_ndts_%i.mat',f0,f1,d,df,phase,ndts), ...
    "continuationarray", "f0", "f1", "d", "df", "phase", "ndts");
end





% save('f0_200_f1_240_d_88_df_1_phase_5.000000e-01_ndts_400.mat')