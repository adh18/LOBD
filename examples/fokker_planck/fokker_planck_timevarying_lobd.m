% Set up the script
clear; close all; 
%cd(fileparts(mfilename('fullpath')))
addpath('../../tensorlab')
addpath('../../functions')

%% Simulate or load the data
if isfile('fokker_planck_data.mat')
    load('fokker_planck_timevarying_data_final.mat')
else
    run('fokker_planck_timevarying_simululations')
end
%% Run the LOBD %1.13950508e-16
tstart = 15; % Time offset in the data
R = 4;       % Number of LOBD bases

%Build the data cell for LOBD 
sols = {sol1(:, tstart:tstart + 800), sol4(:, tstart:tstart + 800)};
usedsols = {'sol1', 'sol2'};
[lobd, output] = LOBD(sols, R, 'useminf', true, 'maxiters', 5000, 'cgiters', 500, 'showevery', 10, 'nonneg', false);

%% Calculate the exact DMD of the same input data
[dmdX, dmdT, omegas, ~, bs] = exactDMD(sols, R, t(2) - t(1));

%% Predictions
newic = sol2(:, tstart);
newsol = sol2(:, tstart:tstart + 800);

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
lobdrelerr = abs(lobdpred - newsol) ./ (1 + abs(newsol));
dmdrelerr = abs(dmdpred - newsol) ./ (1 + abs(newsol));

%%
figure 
subplot(2, 5, 1)
imagesc(t, x, sols{1}); caxis([0, 0.02]); ylabel('Position x')
subplot(2, 5, 2)
imagesc(t, x, sols{2}); caxis([0, 0.02]);
subplot(2, 5, 6);
imagesc(t, x, newsol); caxis([0, 0.02]); xlabel('Time t'); ylabel('Position x')
subplot(2, 5, 7);
imagesc(t, x, lobdpred); caxis([0, 0.02]); xlabel('Time t');
subplot(2, 5, 8);
imagesc(t, x, real(dmdpred)); caxis([0, 0.02]); xlabel('Time t');
subplot(2, 5, 9);
imagesc(t, x, log10(lobdrelerr)); xlabel('Time t'); colormap(gca, 'jet'); caxis([-4, -2]); 
subplot(2, 5, 10);
imagesc(t, x, log10(dmdrelerr)); xlabel('Time t'); colormap(gca, 'jet'); caxis([-4, -2]); 

subplot(2, 5, 3)
plot(x, lobd.factors{1}); ylabel('LOBD spatial basis'); xlabel('Position x')
subplot(2, 5, 4)
plot(x, real(dmdX)); ylabel('DMD spatial basis'); xlabel('Position x')
hold on
plot(x, imag(dmdX), '--'); ylabel('DMD spatial basis'); xlabel('Position x')
subplot(2, 5, 5)
plot(t(1:size(newsol, 2)), lobd.factors{2}); ylabel('LOBD time basis')

%% Save the lobd results
save('fokker_planck_timevarying_lobdresults_final.mat', 'tstart', 'R', 'lobd', 'usedsols', ...
    'lobdpred', 'dmdpred', 'lobdnormerr', 'dmdnormerr', 'lobdrelerr', 'dmdrelerr', ...
    'dmdT', 'dmdX')