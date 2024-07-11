global d
global f
global new_x

data = load('f0_200_f1_240_d_88_df_1_phase_5.000000e-01_ndts_400.mat')

% save("continuationinputantiphase.mat", "d", "f0", "f1", "numf", "new_x")
continuationarray = data.continuationarray;

d = data.d;
f0 = data.f0;
f1 = data.f1;
df = data.df;
numf = (f1-f0)/df;

evalarray = zeros(numf+1, 6);
evecarray = zeros(numf+1, 114, 6);

parfor i=1:numf+1
    f = f0+df*(i-1);
    new_x = continuationarray(:, i);
    fprintf('floquet: starting iteration with f=%i\n',f) ;
    [eval, evec] = LinearStability(f, d, new_x(1), new_x(2:end), 200);
    evalarray(i, :) = eval;
    evecarray(i, :, :) = evec;
end

save(sprintf('linearstability_f0_%i_f1_%i_d_%i_df_%i_phase_0.5.mat',f0,f1,d,df), ...
    "evalarray", "evecarray", "f0", "f1", "d", "df");
