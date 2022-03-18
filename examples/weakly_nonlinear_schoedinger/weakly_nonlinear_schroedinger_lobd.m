% Set up the script
clear; close all; 
cd(fileparts(mfilename('fullpath')))
addpath('../../tensorlab')
addpath('../../functions')

%% Simulate or load the data
if isfile('weakly_nonlinear_schroedinger_data_normalized.mat')
    load('weakly_nonlinear_schroedinger_data.mat')
else
    run('weakly_nonlinear_schroedinger_simulations')
end

%% Run the LOBD
tstart = 1;
R = 7;       % Number of LOBD bases

%Build the data cell for LOBD 
sols = {sol1(:, 1:400), sol2(:, 1:400), sol4(:, 1:400)};
usedsols = {'sol1', 'sol2'};
% For complex data minf must be used and projected coefficients should be
% false (for not : TODO)
[lobd, output, z] = LOBD(sols, R, 'useminf', true, 'maxiters', 20000, 'cgiters', 20000, 'showevery', 10, 'tol', 1e-20, 'projectcoeffs', true);

%%
notdeltas = ~any(abs(lobd.factors{1}) > 0.9);
lobd.factors{1} = lobd.factors{1}(:, notdeltas);
lobd.factors{2} = lobd.factors{2}(:, notdeltas);
%% Calculate the exact DMD of the same input data
R = size(lobd.factors{1}, 2)
[dmdX, dmdT, omega, ~, bs] = exactDMD(sols, R, t(2) - t(1));

%% Predictions
newic = sol3(:, tstart);
newsol = sol3(:, 1:799);
reps = 2;
% Calculate the new coefficients
lobdcfs = lobd.factors{1}' * newic; %./ conj(lobd.factors{2}(1,:)');   % orthogonal projection
dmdcfs = dmdX \ newic;              % least squares fit

% Form the prediction products
lobdpred = LOBDprediction(lobd, lobdcfs, reps);
%lobdpred = cpdgen({lobd.factors{1}, lobd.factors{2}, conj(lobdcfs')});

% Generate longer time DMD matrix: 
dmdTlong = zeros(size(dmdT, 2), size(lobdpred, 2));
tlong = (0:size(lobdpred, 2)-1)*(t(2) - t(1)); % time vector
for iter = 1:size(lobdpred, 2)
    dmdTlong(:,iter) = (exp(omega*tlong(iter)));
end
dmdTlong =dmdTlong';
dmdpred = DMDprediction(dmdX, dmdTlong, dmdcfs);

% Calculate the matrix relative errors
lobdnormerr = norm(lobdpred(101:300, :) - newsol(101:300, :), 'fro')/norm(newsol(101:300, :), 'fro')
dmdnormerr = norm(dmdpred(101:300, :) - newsol(101:300, :), 'fro')/norm(newsol(101:300, :), 'fro')

% Calculate the relative error heatmaps
lobdrelerr = abs(lobdpred - newsol) ./ abs(newsol);
dmdrelerr = abs(dmdpred - newsol) ./ abs(newsol);

%% Save the lobd results
save('weakly_nonlinear_schroedinger_lobdresults.mat', 'tstart', 'R', 'lobd', 'usedsols', ...
    'lobdpred', 'dmdpred', 'lobdnormerr', 'dmdnormerr', 'lobdrelerr', 'dmdrelerr', ...
    'dmdT', 'dmdX')