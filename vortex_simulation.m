%% Simulate BEC vortex dynamics on torus, son!
clear, clc, clf;

% torus parameters;
alpha = 0.7; % alpha = r/R
R = 12; % major radius
r = R*alpha; % minor radius
c = sqrt(R^2 - r^2);

% local curvature
lambda =@(v) c./(R-r*cos(v./r));

% Jacobi theta function parameters
cap = 20; % truncation error
p = exp(-pi*r/c); % nome (for double periodicity)

% complex potential stuff, example intial vortex positions
w1_0 = (-3*pi) + 1i*(-6); % positive vortex
q1 = 1;
w2_0 = (3*pi) + 1i*(0); % negative vortex
q2 = -1;
q = [q1 q2]; % vector of vortex charges
N = length(q); % keeping track of number of vortices

% real and imaginary parts of isothermal coords
u1_0 = real(w1_0); v1_0 = imag(w1_0);
u2_0 = real(w2_0); v2_0 = imag(w2_0);

% integrate these guys
t0 = 0; tf = 1600;
timespan = [t0 tf];
reltolerance = 1e-13; % relative tolerance
abstolerance = 1e-13; % absolute tolerance, for integration

y0 = [u1_0 u2_0 v1_0 v2_0];
options=odeset('Reltol',reltolerance,'Abstol',abstolerance);
[t,y] = ode45('vortex_velocity',timespan,y0,options,N,q,r,R,c,p,cap);

% vortices in isothermal coordinates
U = y(:,1:N); % u-coords
V = y(:,(1+N):2*N); % v-coords
W = U + 1i*V; 

% convert to toroidal and 3D cartesian coords
Phi = U./c;
Theta = 2*atan(sqrt((R+r)/(R-r))*tan(V./(2*r)));
Phi = unwrap(Phi); % reduce mod 2pi
Theta = unwrap(Theta); 

X = (R+r*cos(Theta)).*cos(Phi);
Y = (R+r*cos(Theta)).*sin(Phi);
Z = r*sin(Theta);

%% plot the entire solution 
figure (1);

% toriodal-poloidal
subplot(2,2,1)
plot(Phi./(2*pi),Theta./(2*pi))
grid on
xlabel('$\frac{\phi}{2\pi}$','Interpreter','latex')
ylabel('$\frac{\theta}{2\pi}$','Interpreter','latex')
title('Toroidal-Poloidal','Interpreter','latex')

% 3D Cartesian
subplot(2,2,2)
plot3(X,Y,Z)
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
title('3D Cartesian (physical path)','Interpreter','latex')

% isothermal
subplot(2,2,3)
plot(U,V)
grid on
xlabel('$u = $Re$(w)$','Interpreter','latex')
ylabel('$v = $Im$(w)$','Interpreter','latex')
title('Isothermal','Interpreter','latex')

% 2D Cartesian
subplot(2,2,4)
plot(R*Phi,r*Theta)
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('2D Cartesian','Interpreter','latex')
