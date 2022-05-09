%% Computing fixed points of torus ODE's
% useful after visualizing bits of phase space
clear, clc;
format long
figure(1)

% torus parameters;
alpha0 = 0.4; % alpha = r/R
R = 12; % major radius
r = R*alpha0; % minor radius
c = sqrt(R^2 - r^2);

% Jacobi theta function parameters
cap = 12; % truncation error
p = exp(-pi*r/c); % nome

% intial (educated) guess of fixed points positions
u1 = c*(asin(3.9874/16.8)); v1 = 0;   % positive vortex
u2 = c*(asin(-3.9874/16.8)+2*pi); v2 = 0; % negative vortex
q1 = 1; q2 = -1;
q = [q1 q2]; % vortex charges
N = length(q); % number of vortices

%% Solve nonlinear least squares problem to estimate fixed points
% Levenberg-Marquard fixed point iteration
% z = [u1, u2, v1, v2]

z0 = [u1, u2, v1, v2];   % z0 is informed guess of FP
options = optimoptions('lsqnonlin','Display','iter');
options.Algorithm = 'levenberg-marquardt';
[zstar,resnorm,residual,exitflag,output] =...
    lsqnonlin(@skroot,z0,[],[],options);
fprintf('zstar is %0.8f\n',zstar)
disp(output);

%% S(z) function
function S = skroot(z)
% Formulates torus RHS = 0 as a nonlinear lsq problem
% takes in:
% z = [u1, u2, v1, v2]
% params = [N,q,r,R,c,p,cap]
N = 2; q1 = 1; q2 = -1; q = [q1 q2]; R = 12;alpha0 = 0.4; r = R*alpha0;
c = sqrt(R^2 - r^2); p = exp(-pi*r/c); cap = 12;
% spits out:
% S = (dydt)^(1/2)
    dydt = vortex_velocity([],z,[],N,q,r,R,c,p,cap);
    S = sqrt(dydt);
end






