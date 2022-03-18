% Set up the script
clear; close all; 
%cd(fileparts(mfilename('fullpath')))
addpath('../../tensorlab')
addpath('../../functions')

%%

% Simulation parameters
Nx = 500;                   %Number of spatial points
Nt = 500;                   %Number of time points
t = linspace(0,20,Nt);      %t range
x = linspace(-10,10,Nx);    %x range

ic1 = @(x) (5*x.^4 - x + 2*x.^3 + 0.4*x.^2).*exp(-0.5*x.^2) / 17.6964;
ic3 = @(x) (2*x.^3 + 6*x.^2 + 1.5*x + 5 + x.^4).*exp(-0.5*x.^2) / 15.2778;
ic2 = @(x) (3*x.^4 + x - sin(5)*x.^3 + pi*x.^2 - 4).*exp(-0.5*x.^2) / 12.6934;
%ic3 = @(x) (x.^3 - 4*x.^2 + (1/3)*x - 3 + rand()*x.^4).*exp(-0.5*x.^2) / 7.7161;
%ic5 = sin(4*x - 40).*exp(-0.5*x.^2)/norm(sin(4*x - 40).*exp(-0.5*x.^2));
%ic6 = cos(x).^2.*exp(-0.5*x.^2);

%%
% Simulate the PDEs
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
%opts = odeset();
sol1 = pdepe(0,@(x,t,u,dudx) qhofun(x,t,u,dudx),ic1,@bcfun,x,t,opts)';
disp('Done 1')
sol2 = pdepe(0,@(x,t,u,dudx) qhofun(x,t,u,dudx),ic2,@bcfun,x,t,opts)';
disp('Done 2')
sol3 = pdepe(0,@(x,t,u,dudx) qhofun(x,t,u,dudx),ic3,@bcfun,x,t,opts)';
disp('Done 3')
%sol4 = pdepe(0,@(x,t,u,dudx) fpefun(x,t,u,dudx,D),ic4,@bcfun,x,t,opts)';

%%
% Save the data 
save('quantum_harmonic_oscillator_data.mat', 'sol1', 'sol2', 'sol3', 'x', 't', 'Nt', 'Nx');

%% Local functions for the FPE
function [c,f,s] = qhofun(x,t,u,dudx)
c = 1;
f = 0.5*1i*dudx;
s = -0.5*1i*x.^2.*u;
end

function [pL,qL,pR,qR] = bcfun(xL,uL,xR,uR,t)
pL = uL;
qL = 0;
pR = uR;
qR = 0;
end