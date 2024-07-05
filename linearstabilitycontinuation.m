global d
global f
global new_x

data = load('f0_100_f1_200_d_88_df_1_phase_0.5.mat')

% save("continuationinputantiphase.mat", "d", "f0", "f1", "numf", "new_x")
continuationarray = data.continuationarray;

d = 88;
f0 = 100;
f1 = 200;
numf = 100;

df = (f1-f0)/numf;

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

save(sprintf('linearstability_f0_%i_f1_%i_d_%i_df_%i_phase_0.5',f0,f1,d,df), ...
    "evalarray", "evecarray", "f0", "f1", "d", "df");
