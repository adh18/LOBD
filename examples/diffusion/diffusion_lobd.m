% Set up the script
clear; close all; 
%cd(fileparts(mfilename('fullpath')))
addpath('../../tensorlab')
addpath('../../functions')

%% Simulate or load the data
if isfile('diffusion_data.mat')
    load('diffusion_data.mat')
else
    run('diffusion_simululations')
end

%% Run the LOBD
tstart = 40; % Time offset in the data
R = 6;      % Number of LOBD bases

%tstart = 40; % Time offset in the data
%R = 6;      % Number of LOBD bases

%Build the data cell for LOBD 
sols = {sol1(:, tstart:end), sol2(:, tstart:end)};
usedsols = {'sol1', 'sol2'};
rank([sols{1} sols{2}])
[lobd, output, z1] = LOBD(sols, R, 'useminf', false, 'maxiters', 200, 'cgiters', 5000, 'showevery', 1);
lobd1 = lobd;

%%
lobd = lobd1;

%% Calculate the exact DMD of the same input data
[dmdX, dmdT, ~, ~, bs] = exactDMD(sols, R, t(2) - t(1));

%% Predictions
newic = sol3(:, tstart);
newsol = sol3(:, tstart:end);

% Calculate the new coefficients
lobdcfs = lobd.factors{1}'*newic;   % orthogonal projection
dmdcfs = dmdX \ newic;              % least squares fit

% Form the prediction products
lobdpred = LOBDprediction(lobd, lobdcfs);
dmdpred = DMDprediction(dmdX, dmdT, dmdcfs);

% Calculate the matrix relative errors
lobdnormerr = norm(lobdpred - newsol, 'fro')/norm(newsol, 'fro')
dmdnormerr = norm(dmdpred - newsol, 'fro')/norm(newsol, 'fro')

% Calculate the relative error heatmaps
lobdrelerr = abs(lobdpred - newsol) ./ abs(newsol);
dmdrelerr = abs(dmdpred - newsol) ./ abs(newsol);

%% Save the lobd results
save('diffusion_lobdresults.mat', 'tstart', 'R', 'lobd', 'usedsols', ...
    'lobdpred', 'dmdpred', 'lobdnormerr', 'dmdnormerr', 'lobdrelerr', 'dmdrelerr', ...
    'dmdT', 'dmdX')