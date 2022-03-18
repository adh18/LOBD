% Set up the script
clear; close all; 
cd(fileparts(mfilename('fullpath')))
addpath('../../tensorlab')
addpath('../../functions')

%% Generate the test data
% We simulate the forced diffusion equation using pdepe. The definitions
% are at the end of the script in diffsionfun and bcfun 

% Simulation parameters
Nx = 200;                   %Number of spatial points
Nt = 300;                   %Number of time points
t = linspace(0, 1, Nt + 1); %t range
x = linspace(0, 1, Nx);     %x range
D = 0.1;                    %Diffusion constant

% Functions for the initial conditions
ic1 = @(xfun) exp(-50*(xfun - 0.4).^2) + exp(-50*(xfun + 0.4).^2);
ic2 = @(xfun) 0.25*exp(-500*(xfun - 0.6).^2) + 0.75*exp(-100*(xfun - 0.3).^2);
ic3 = @(xfun) 1./cosh(40*(xfun - 0.4));

% Simulate the PDEs
opts = odeset('RelTol',2.5e-14, 'AbsTol',1e-20); %PDE simuation parameters
sol1 = pdepe(0,@(x,t,u,dudx) diffusionfun(x,t,u,dudx,D),ic1,@bcfun,x,t,opts)';
sol2 = pdepe(0,@(x,t,u,dudx) diffusionfun(x,t,u,dudx,D),ic2,@bcfun,x,t,opts)';
sol3 = pdepe(0,@(x,t,u,dudx) diffusionfun(x,t,u,dudx,D),ic3,@bcfun,x,t,opts)';

%%
% Save the data 
save('diffusion_data.mat', 'sol1', 'sol2', 'sol3', 'x', 't', 'Nt', 'Nx');

% Plot a check of the simulation
fig = figure('units', 'normalized', 'position', [0.1, 0.1, 0.6, 0.2]);
clim = [0, 0.3]; % Colormap limits

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

%% Local functions for simulating the diffusion equation
function [c,f,s] = diffusionfun(x,t,u,dudx,D)
c = 1;
f = D*dudx;
s = 0;
end

% Enforce Dirichlet BCs
function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t)
pL = uL;
qL = 0;
pR = uR;
qR = 0;
end