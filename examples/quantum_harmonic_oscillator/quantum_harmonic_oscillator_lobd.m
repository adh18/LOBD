% Set up the script
clear; close all; 
%cd(fileparts(mfilename('fullpath')))
addpath('../../tensorlab')
addpath('../../functions')

%% Simulate or load the data
if isfile('quantum_harmonic_oscillator_data_normalized.mat')
    load('quantum_harmonic_oscillator_data.mat')
else
    run('quantum_harmonic_oscillator_simulations')
end

%% Run the LOBD
tstart = 1;
R = 7;       % Number of LOBD bases

%Build the data cell for LOBD 
sols = {sol1, sol2};
usedsols = {'sol1', 'sol2'};
% For complex data minf must be used and projected coefficients should be
% false (for not : TODO)
[lobd, output, z] = LOBD(sols, R, 'useminf', true, 'maxiters', 100000, 'cgiters', 5000, 'showevery', 10, 'tol', 1e-30, 'projectcoeffs', true);

%%
notdeltas = ~any(abs(lobd.factors{1}) > 0.9);
lobd.factors{1} = lobd.factors{1}(:, notdeltas);
lobd.factors{2} = lobd.factors{2}(:, notdeltas);
%% Calculate the exact DMD of the same input data
R = size(lobd.factors{1}, 2)
[dmdX, dmdT, ~, ~, bs] = exactDMD(sols, R, t(2) - t(1));

%% Predictions
newic = sol3(:, tstart);
newsol = sol3(:, tstart:end);

% Calculate the new coefficients
lobdcfs = lobd.factors{1}'*newic ./ conj(lobd.factors{2}(1,:)');   % orthogonal projection
dmdcfs = dmdX \ newic;              % least squares fit

% Form the prediction products
lobdpred = LOBDprediction(lobd, lobdcfs);
dmdpred = DMDprediction(dmdX, dmdT, dmdcfs);

% Calculate the matrix relative errors
lobdnormerr = norm(lobdpred(151:350, :) - newsol(151:350, :), 'fro')/norm(newsol(151:350, :), 'fro')
dmdnormerr = norm(dmdpred(151:350, :) - newsol(151:350, :), 'fro')/norm(newsol(151:350, :), 'fro')

% Calculate the relative error heatmaps
lobdrelerr = abs(lobdpred - newsol) ./ abs(newsol);
dmdrelerr = abs(dmdpred - newsol) ./ abs(newsol);

%% Save the lobd results
save('quantum_harmonic_oscillator_lobdresults.mat', 'tstart', 'R', 'lobd', 'usedsols', ...
    'lobdpred', 'dmdpred', 'lobdnormerr', 'dmdnormerr', 'lobdrelerr', 'dmdrelerr', ...
    'dmdT', 'dmdX')