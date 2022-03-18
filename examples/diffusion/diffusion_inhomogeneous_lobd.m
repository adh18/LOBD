% Set up the script
clear; close all; 
cd(fileparts(mfilename('fullpath')))
addpath('../../tensorlab')
addpath('../../functions')

%% Simulate or load the data
if isfile('diffusion_inhomogeneous_data.mat')
    load('diffusion_inhomogeneous_data.mat')
else
    run('diffusion_inhomogeneous_simululations')
end

%% Run the LOBD
tstart = 7; % Time offset in the data
R = 6;      % Number of LOBD bases

%Build the data cell for LOBD 
sols = {sol3(:, tstart:end), sol2(:, tstart:end)};
usedsols = {'sol3', 'sol2'};
[lobd, output] = LOBD(sols, R, 'useminf', false, 'maxiters', 1000, 'cgiters', 5000);

%% Calculate the exact DMD of the same input data
[dmdX, dmdT, ~, ~, bs] = exactDMD(sols, R, t(2) - t(1));

%% Predictions
newic = sol1(:, tstart);
newsol = sol1(:, tstart:end);

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
save('diffusion_inhomogeneous_lobdresults.mat', 'tstart', 'R', 'lobd', 'usedsols', ...
    'lobdpred', 'dmdpred', 'lobdnormerr', 'dmdnormerr', 'lobdrelerr', 'dmdrelerr', ...
    'dmdT', 'dmdX')