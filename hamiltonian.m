function [energy,classic,curve,quantum] = hamiltonian(U,V,N,q,p,c,r,R,cap)
% This function computes the total energy of a system of vortices, where:

% W = U + iV = wortex positions vector (in isothermal coorindates)
% N = number of vortices

% q = vector of vortex charges (should sum to 0)
% r,R,c = torus parameters
% p = nome of first Jacobi theta function (should satisfy |p| < 1)
% cap = truncation term of sum in first Jacobi theta function

% saving on argin's 
W = U + 1i*V; 

% local curvature
L =@(v) c./(R-r*cos(v./r));

% initialize sum variables
classic = 0; curve = 0; quantum = 0;

% loop over each vortex pair
for n = 1:N
    % curvature contribution
    curve = curve + 2*pi*log(L(V(:,n)));

    for m = 1:N
        % quantum contribution
        quantum = quantum + (1/(r*c)).*q(n).*q(m).*U(:,n).*U(:,m);

        if (m ~= n)
            % classical contribution
            thetus = (W(:,n) - W(:,m))./(2*c);
            classic = classic -...
                2*pi*(q(n).*q(m).*log(abs(jacobitheta1(thetus,p,cap))));
        end
    end
end
energy = classic + curve + quantum;
end