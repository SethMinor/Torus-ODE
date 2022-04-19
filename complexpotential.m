function Fofw = complexpotential(w,w1,w2,r,c,p,cap)
% This function computes the complex potential (evaluated at w) given 
% by Fetter in eq (14), where:
% w = u + vi = complex number to evaluate potential at
% w1 = u1 + v1i = position of first (positive) vortex (isothermal)
% w2 = u2 + v2i = position of second (negative) vortex (isothermal)
% r and c = torus parameters
% p = nome of first Jacobi theta function (should satisfy |p| < 1)
% cap = truncation error of sum in first Jacobi theta function
% Fofw = value of the complex potential at complex coord w = u + vi

topinput = (w-w1)/(2*c);
utop = real(topinput);
vtop = imag(topinput);

botinput = (w-w2)/(2*c);
ubot = real(botinput);
vbot= imag(botinput);

theta1 = jacobitheta1(utop,vtop,p,cap);
theta2 = jacobitheta1(ubot,vbot,p,cap);

linearterm = ((real(w1-w2))./(2*pi*r*c))*w;

Fofw = (log(theta1)-log(theta2)) - linearterm;
end