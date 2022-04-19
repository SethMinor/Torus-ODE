function dydt = vortex_velocity(~,y0,~,N,q,r,R,c,p,cap)
%% This function computes the physical velocity of a set of vortices,
% using the relation given by eq (35) in Fetter, where:

% y0 = initial conditions vector, e.g. y0 = [u1_0 u2_0 v1_0 v2_0]
%       (in isothermal coorindates)
% N = number of vortices

% q = vector of vortex charges (should sum to 0)
% r,R,c = torus parameters
% p = nome of first Jacobi theta function (should satisfy |p| < 1)
% cap = truncation term of sum in first Jacobi theta function

% omega = Vy + iVx = complex physical velocity of vortex cores
% dydt = right-hand side of ODE

%% right-hand side of ODE
% real and imaginary parts of isothermal vortex coords
y0 = y0(:);
U = y0(1:N);
V = y0((N+1):2*N);
W = U+1i*V;

% curvature terms
L = c./(R-r*(cos(V./r)));
Lprime = -c*sin(V./r)./((R-r*cos(V./r)).^2);

% initialize omega vector
omega = zeros(N,1);

% loop over each vortex
for n=1:N
    % add curvature and quantum contributions
    curve = (1i/2)*q(n)*(Lprime(n)./L(n));
    quantum = -(1/(2*pi*r*c))*(q(n).*U(n));
    omega(n) = omega(n) + curve + quantum;

    % add induced velocity contributions (from other vortices)
    for m=1:N
        if (m ~= n)
            T = jacobitheta1((W(n)-W(m))./(2*c),p,cap); % jacobi theta
            Tprime = Djacobitheta1((W(n)-W(m))./(2*c),p,cap);
            fmn = 1/(2*c)*(Tprime./T)-(1/(2*pi*r*c))*U(m);
            omega(n) = omega(n) + q(m).*fmn;
        end
    end

    omega(n) = (1./L(n)).*omega(n);
end

vdots = real(omega)./L;
udots = imag(omega)./L;

dydt = [udots; vdots];
end