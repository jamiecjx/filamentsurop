global d
global f

data = load("continuationinput.mat")

% d = 2.2*N;
% f0 = 100;
% f1 = 200;
% numf = 100;
% new_x = xt0;

d = data.d;
f0 = data.f0;
f1 = data.f1;
numf = data.numf;
new_x = data.xt0;

df = (f1-f0)/numf;

continuationarray = zeros(size(new_x, 1), numf+1);
continuationarray(:,1) = new_x;

for i=1:numf
    f = f0+df*i;
    fprintf('continuation: starting iteration with f=%i',f) ;
    JFNK
    continuationarray(:,i+1) = new_x;
end

save(sprintf('f0_%i_f1_%i_d_%i_df_%i',f0,f1,d,df), "continuationarray");
