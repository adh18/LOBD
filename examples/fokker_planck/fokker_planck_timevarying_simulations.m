% Set up the script
clear; close all; 
%cd(fileparts(mfilename('fullpath')))
addpath('../../tensorlab')
addpath('../../functions')

%% Generate the test data
% We simulate the forced diffusion equation using pdepe. The definitions
% are at the end of the script in diffsionfun and bcfun 

% Simulation parameters
Nx = 300;                   %Number of spatial points
Nt = 1000;                  %Number of time points
t = linspace(0,100,Nt);     %t range
x = linspace(-3,3,Nx);      %x range
D = 0.45;                    %Diffusion constant

% Functions for the initial conditions
ic1 = @(xfun) (exp(-10*(xfun - 0.8).^2) + exp(-10*(xfun + 1.1).^2))/69.9223;
ic2 = @(xfun) 5*exp(-10*(xfun - 0.1).^2)/ 174.8057; %+ 0.75*exp(-7*(xfun + 1.4).^2);
ic3 = @(xfun) (5*exp(-20*(xfun - 2).^2) + 5*exp(-20*(xfun + 2).^2) + 5*exp(-20*(xfun).^2))/524.4159/0.7071;
ic4 = @(xfun) 5*exp(-20*(xfun - 0.4).^2)/ 174.8057/ 0.7071;

% Simulate the PDEs
opts = odeset('RelTol',2.25e-14,'AbsTol',1e-18);
sol1 = pdepe(0,@(x,t,u,dudx) fpefun(x,t,u,dudx,D),ic1,@bcfun,x,t,opts)';
disp('Done 1')
sol2 = pdepe(0,@(x,t,u,dudx) fpefun(x,t,u,dudx,D),ic2,@bcfun,x,t,opts)';
disp('Done 2')
sol3 = pdepe(0,@(x,t,u,dudx) fpefun(x,t,u,dudx,D),ic3,@bcfun,x,t,opts)';
disp('Done 3')
sol4 = pdepe(0,@(x,t,u,dudx) fpefun(x,t,u,dudx,D),ic4,@bcfun,x,t,opts)';
disp('Done 4')
%%
% Save the data 
save('fokker_planck_timevarying_data_final.mat', 'sol1', 'sol2', 'sol3', 'sol4', 'x', 't', 'Nt', 'Nx');

%%
% Plot a check of the simulation
fig = figure('units', 'normalized', 'position', [0.1, 0.1, 0.6, 0.2]);
clim = [0, 0.015]; % Colormap limits

% Plot the initial conditions
subplot(1, 4, 1)
plot(x, ic1(x), 'LineWidth', 1)
hold on; title('Initial Conditions'); xlabel('Position x'); ylabel('Density');
plot(x, ic2(x), 'LineWidth', 1)
plot(x, ic3(x), 'LineWidth', 1)
leg = legend('Sample 1', 'Sample 2', 'Sample 3', 'Box', 'off'); leg.ItemTokenSize = [10, 18];

% Plot the samples
subplot(1, 4, 2)
imagesc(t, x, sol1); ax = gca; title('Sample 1'); xlabel('Time t'); caxis(clim); ylabel('Position x')
subplot(1, 4, 3)
imagesc(t, x, sol2); ax = gca; title('Sample 2'); xlabel('Time t'); caxis(clim); ax.YTickLabels = []; 
subplot(1, 4, 4)
imagesc(t, x, sol3); ax = gca; title('Sample 3'); xlabel('Time t'); caxis(clim); ax.YTickLabels = [];
%saveas(fig, 'diffusion_homogeneous_simulation_data', 'pdf')

%% Local functions for the FPE
function [c,f,s] = fpefun(x,t,u,dudx,D)
c = 1;
f = D  * (1 + 0.2*sin(0.2*t))^2*dudx;
s = (3*x.^2 - 3*1).*u + (x.^3 - 3*x).*dudx;
end

function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t)
pL = uL;
qL = 0;
pR = uR;
qR = 0;
end