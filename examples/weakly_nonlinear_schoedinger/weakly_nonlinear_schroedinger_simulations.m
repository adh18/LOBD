% Set up the script
clear; close all; 
%cd(fileparts(mfilename('fullpath')))
addpath('../../tensorlab')
addpath('../../functions')

%%

% Simulation parameters
Nx = 400;                   %Number of spatial points
Nt = 800;                   %Number of time points
t = linspace(0,2,Nt);      %t range
x = linspace(-10,10,Nx);    %x range

xfun = x;
%ic1 = @(xfun) (5*xfun.^4 - xfun + 2*xfun.^3 + 0.4*xfun.^2).*exp(-0.5*xfun.^2) / 17.6964;
%ic2 = @(xfun) (xfun.^3 + 3*xfun.^2 + 0.5*xfun + 1 - xfun.^4).*exp(-0.5*xfun.^2) / 3.6001;
%ic3 = @(xfun) (3*xfun.^4 + xfun - sin(5)*xfun.^3 + pi*xfun.^2 - 4).*exp(-0.5*xfun.^2) / 12.6934;
ic4 = @(xfun) (-0.3*xfun.^3 - 1.5*xfun.^2 + 0.2*xfun - 0.3 + xfun.^4).*exp(-0.5*xfun.^2) / 2.2332;
ic1 = @(x) (5*x.^4 - x + 2*x.^3 + 0.4*x.^2 - 1).*exp(-0.5*x.^2) / 17.3474;
ic3 = @(x) (2*x.^3 + 6*x.^2 + 1.5*x + 5 + x.^4).*exp(-0.5*x.^2) / 15.2778;
ic2 = @(x) (3*x.^4 + x - sin(5)*x.^3 + pi*x.^2 - 4).*exp(-0.5*x.^2) / 12.6934;

%%
% Simulate the PDEs
opts = odeset('RelTol',2.5e-13,'AbsTol',1e-14);
sol1 = pdepe(0,@(x,t,u,dudx) qhofun(x,t,u,dudx),ic1,@bcfun,x,t,opts)';
1
sol2 = pdepe(0,@(x,t,u,dudx) qhofun(x,t,u,dudx),ic2,@bcfun,x,t,opts)';
2
sol3 = pdepe(0,@(x,t,u,dudx) qhofun(x,t,u,dudx),ic3,@bcfun,x,t,opts)';
3
sol4 = pdepe(0,@(x,t,u,dudx) qhofun(x,t,u,dudx),ic4,@bcfun,x,t,opts)';
4

%%
% Save the data 
save('weakly_nonlinear_schroedinger_data.mat', 'sol1', 'sol2', 'sol3', 'sol4', 'x', 't', 'Nt', 'Nx');

%% Local functions for the NLS
function [c,f,s] = qhofun(x,t,u,dudx)
c = 2*i;
f = -dudx;
s = -0.5 * u.*abs(u).^2 + x.^2.*u;
end

function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t)
pL = uL;
qL = 0;
pR = uR;
qR = 0;
end

