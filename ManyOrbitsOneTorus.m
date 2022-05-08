%% Visualize bits of torus ODE phase space
% plot several orbits on the same torus
clear, clc, clf;
figure (1);

% torus parameters;
alpha0 = 0.4; % alpha = r/R
R = 12; % major radius
r = R*alpha0; % minor radius
c = sqrt(R^2 - r^2);

% Jacobi theta function parameters
cap = 12; % truncation error
p = exp(-pi*r/c); % nome

% put the torus in the background
utorus = linspace(0,2*pi);
vtorus = linspace(0,2*pi);
[Utorus,Vtorus] = meshgrid(utorus,vtorus);
Xtorus = (R+r.*cos(Vtorus)).*cos(Utorus);
Ytorus = (R+r.*cos(Vtorus)).*sin(Utorus);
Ztorus = r.*sin(Vtorus);
surf(Xtorus, Ytorus, Ztorus,'FaceAlpha',0.6);
colormap('bone'); shading interp;

% plot settings
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
title("Horizontal dipoles $(\alpha=$ "+alpha0+")",...
    'interpreter','latex','FontSize',12)
axis equal

% some pretty colors (for plotting)
bluey = [0 0.4470 0.7410];
orangu = [0.8500 0.3250 0.0980];

% integration settings
t0 = 0; tf = 500;
timespan = [t0 tf];
reltolerance = 1e-12; % relative tolerance
abstolerance = 1e-12; % absolute tolerance 

%% fixing some parameters and varying others while plotting

% define range of varying parameter
num_plots = 160; % number of trajectories on torus
sneb = linspace(0, 53.8, num_plots);

hold on
for j = 1:num_plots
    disp('j is ')
    disp(j)
    greb = sneb(j);

    u1 = greb; v1 = 0; % positive vortex
    u2 = -greb; v2 = 0; % negative vortex
    q1 = 1; q2 = -1;
    q = [q1 q2]; % vortex charges
    N = length(q); % number of vortices
    
    y0 = [u1 u2 v1 v2];
    disp('y0 is ')
    disp(y0)
    options=odeset('Reltol',reltolerance,'Abstol',abstolerance);
    [t,y] = ode45('vortex_velocity',timespan,y0,options,N,q,r,R,c,p,cap);
    disp('y is ')
    disp(y)
    
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
    
    plot3(X(:,1),Y(:,1),Z(:,1),'LineWidth',0.5,'Color',orangu)
    plot3(X(:,2),Y(:,2),Z(:,2),'LineWidth',0.5,'Color',bluey)
end
hold off

disp('done')




