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
cap = 12; % truncation error
p = exp(-pi*r/c); % nome (for double periodicity)

% complex potential stuff, example intial vortex positions
w1_0 = (-2.3*pi) + 1i*(-3); % positive vortex
q1 = 1;
w2_0 = (1*pi) + 1i*(1); % negative vortex
q2 = -1;
q = [q1 q2]; % vector of vortex charges
N = length(q); % keeping track of number of vortices
% create vortex class in matlab?

% real and imaginary parts of isothermal coords
u1_0 = real(w1_0); v1_0 = imag(w1_0);
u2_0 = real(w2_0); v2_0 = imag(w2_0);

% integrate these guys
t0 = 0; tf = 1600;
timespan = [t0 tf];
reltolerance = 1e-12; % relative tolerance
abstolerance = 1e-12; % absolute tolerance, for integration
% this also controls precision of total energy

y0 = [u1_0 u2_0 v1_0 v2_0];
options=odeset('Reltol',reltolerance,'Abstol',abstolerance);
[t,y] = ode45('vortex_velocity',timespan,y0,options,N,q,r,R,c,p,cap);

% Numerical Jacobian, Eigenvalues
NN=length(y0);
Jnum=[];
if(exist('dh')==0) dh=1e-4; end

for jj=1:NN
 pert=zeros(NN,1);
 pert(jj)=dh*norm(y);
 ypert=y0(:)+pert;
 rhs = vortex_velocity([],ypert,[],N,q,r,R,c,p,cap);
 drhs = rhs./(dh*norm(y));
 Jnum(:,jj) = drhs;
end
disp(Jnum);
[XX,DD] = eig(Jnum);

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
hold on
plotwrapped(Phi(:,1),Theta(:,1),1, [-pi pi],[-pi pi], 0.05,bluey)
plotwrapped(Phi(:,2),Theta(:,2),1, [-pi pi],[-pi pi], 0.05,orangu)
hold off
xlabel('$\phi$','Interpreter','latex')
ylabel('$\theta$','Interpreter','latex')
title('Toroidal-Poloidal','Interpreter','latex','FontSize',12)
xlim([-pi pi])

% 3D Cartesian
subplot(2,2,2)
% plot the torus
% (see Sathyanarayan Rao,
% 'Visualizing a Toroidal Surface (Torus) in Matlab')
utorus = linspace(0,2*pi);
vtorus = linspace(0,2*pi);
[Utorus,Vtorus] = meshgrid(utorus,vtorus);
Xtorus = (R+r.*cos(Vtorus)).*cos(Utorus);
Ytorus = (R+r.*cos(Vtorus)).*sin(Utorus);
Ztorus = r.*sin(Vtorus);
surf(Xtorus, Ytorus, Ztorus);
colormap('gray')
shading interp;
hold on
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

%% compute and plot Hamiltonian
[energy,classic,curve,quantum] = hamiltonian(U,V,N,q,p,c,r,R,cap);

figure (3);

energy_time = linspace(t0,tf,length(energy));

% plot relative difference (E(t)-E(0))/E(0)
subplot(1,2,1)
plot(energy_time,((energy-energy(1))./energy(1)))
title('Energy of solution','Interpreter','latex')
grid on
xlabel('$t$','Interpreter','latex')
ylabel('$\frac{H(t)-H_0}{H_0}$','Interpreter','latex')

% plot each energy contribution
subplot(1,2,2)
plot(energy_time,energy,'-',energy_time,classic,'--')
hold on
plot(energy_time,curve,'--',energy_time,quantum,'--')
grid on
xlabel('$t$','Interpreter','latex')
ylabel('Energie coontributions','Interpreter','latex')
legend('Total, $H$','Classic','Curvature',...
    'Quantum','Interpreter','latex')
    
    
%% Movie maker Code - Use Independently from Above
clear, clc, clf;

% torus parameters;
alpha0 = 0.4; % alpha = r/R
R = 12; % major radius
r = R*alpha0; % minor radius
c = sqrt(R^2 - r^2);

% Jacobi theta function parameters
cap = 12; % truncation error
p = exp(-pi*r/c); % nome (for double periodicity)

% complex potential stuff, example intial vortex positions
% (enter vortex positions, or enable random seed)

% random seed (enable rng('default') to recover results)
%rng('default');
% randoU = (2*pi*c).*rand(1,2); %- pi*c;
% randoV = (2*pi*r).*rand(1,2); %- pi*r;
% 
% % enter vortex positions (real and imaginary parts, in isothermal coords)
% u1_0 = randoU(1); v1_0 = randoV(1);
% u2_0 = randoU(2); v2_0 = randoV(2);

u1_0 = 44.8562; v1_0 = 19.5356;
u2_0 = 25.5647; v2_0 = 13.5995;

% w1 = (u1 + iv1) is the positive vortex
% w2 = (u2 + iv2) is the negative vortex
q1 = 1; q2 = -1;
q = [q1 q2]; % vector of vortex charges
N = length(q); % keeping track of number of vortices

% integrate these guys
t0 = 0; tf = 1500;
timespan = [t0 tf];
reltolerance = 1e-12; % relative tolerance (1e-12)
abstolerance = 1e-12; % absolute tolerance, for integration (1e-12)
% this also controls precision of total energy

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

% plot the entire solution 
% some color stuff
bluey = [0 0.4470 0.7410];
orangu = [0.8500 0.3250 0.0980];



phimod = wrapToPi(Phi);
thetamod = wrapToPi(Theta);
tol = 0.15;
M = length(phimod);

% toriodal-poloidal
% * also works ---> plot(Phi./(2*pi),Theta./(2*pi),'LineWidth',1)
subplot(2,2,1)
%subplot(2,1,1)

% v = VideoWriter('TorusMovie2.avi');
% v.Quality = 95;
% open(v);
hold on
for m = 1:M
    if (m == 1) % plot starting position
        if (abs(phimod(m,1)-phimod(m+1,1)) < (2*pi-tol)) &&...
                (abs(thetamod(m,1)-thetamod(m+1,1)) < (2*pi-tol)) &&...
                    (abs(phimod(m,2)-phimod(m+1,2)) < (2*pi-tol)) &&...
                        (abs(thetamod(m,2)-thetamod(m+1,2)) < (2*pi-tol))
            % in this case, the plot condition is satisfied
            plot(phimod(m:m+1,1),thetamod(m:m+1,1),...
                'LineWidth',1,'Color',bluey)
            plot(phimod(m:m+1,2),thetamod(m:m+1,2),...
                'LineWidth',1,'Color',orangu)
                grid on
                xlim([-pi pi])
                ylim([-pi pi])
                xlabel('$\phi$','Interpreter','latex')
                ylabel('$\theta$','Interpreter','latex')
                title('Toroidal-Poloidal','Interpreter','latex','FontSize',12)
        end
        plot(phimod(m,1),thetamod(m,1),'ok','MarkerSize',8)
        plot(phimod(m,2),thetamod(m,2),'ok','MarkerSize',8)
            grid on
            xlim([-pi pi])
            ylim([-pi pi])
            xlabel('$\phi$','Interpreter','latex')
            ylabel('$\theta$','Interpreter','latex')
            title('Toroidal-Poloidal','Interpreter','latex','FontSize',12)
    elseif (m > 1) && (m < M) % do regular plotting
        if (abs(phimod(m,1)-phimod(m+1,1)) < (2*pi-tol)) &&...
                (abs(thetamod(m,1)-thetamod(m+1,1)) < (2*pi-tol)) &&...
                    (abs(phimod(m,2)-phimod(m+1,2)) < (2*pi-tol)) &&...
                        (abs(thetamod(m,2)-thetamod(m+1,2)) < (2*pi-tol))
            % in this case, the plot condition is satisfied
            plot(phimod(m:m+1,1),thetamod(m:m+1,1),...,
                'LineWidth',1,'Color',bluey)
            plot(phimod(m:m+1,2),thetamod(m:m+1,2),...,
                'LineWidth',1,'Color',orangu)
                grid on
                xlim([-pi pi])
                ylim([-pi pi])
                xlabel('$\phi$','Interpreter','latex')
                ylabel('$\theta$','Interpreter','latex')
                title('Toroidal-Poloidal','Interpreter','latex','FontSize',12)
        end
    else % plot finishing position
        plot(phimod(m,1),thetamod(m,1),...
            'ok','MarkerFaceColor','k','MarkerSize',8)
        plot(phimod(m,2),thetamod(m,2),...
            'ok','MarkerFaceColor','k','MarkerSize',8)
                grid on
                xlim([-pi pi])
                ylim([-pi pi])
                xlabel('$\phi$','Interpreter','latex')
                ylabel('$\theta$','Interpreter','latex')
                title('Toroidal-Poloidal','Interpreter','latex','FontSize',12)
    end
% pause(0.000001)
% set(gcf);
% frame = getframe(gcf);
% writeVideo(v,frame);
end

%close(v);

% hold off
% grid on
% xlim([-pi pi])
% ylim([-pi pi])

% 3D Cartesian
subplot(2,2,[3 4])
%subplot(2,1,2)
% plot the torus
% (see Sathyanarayan Rao: 'Visualizing a Torus in Matlab')
utorus = linspace(0,2*pi);
vtorus = linspace(0,2*pi);
[Utorus,Vtorus] = meshgrid(utorus,vtorus);
Xtorus = (R+r.*cos(Vtorus)).*cos(Utorus);
Ytorus = (R+r.*cos(Vtorus)).*sin(Utorus);
Ztorus = r.*sin(Vtorus);
surf(Xtorus, Ytorus, Ztorus,'FaceAlpha',0.3);
colormap('bone'); shading interp;
hold on
plot3(X,Y,Z,'LineWidth',1.5)
hold off
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
title('3D Cartesian (physical path)','Interpreter','latex','FontSize',12)
axis equal

% isothermal
subplot(2,2,2)
plot(U./(2*pi*c),V./(2*pi*r),'LineWidth',0.8)
grid on
xlabel('$u = $Re$(w)$','Interpreter','latex')
ylabel('$v = $Im$(w)$','Interpreter','latex')
title('Isothermal','Interpreter','latex','FontSize',12)
