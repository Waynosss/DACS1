clear
close
%% Variables
M = 15e6;
V = 1.5e6;
R = 3;

Ex = 69e9; 
Ey = 69e9;
Gxy =26;
vxy = 0.29;

Xt = 410e6;
Xc = 430e6;
Yt = 4006;
Yc = 430e6;
S = 230e6;

rho = 2770;

thick = 7.6e-3;

%% Calcs
bucklingcritmoment = sigmacrit(Ex, vxy, R, thick) * pi * R^2 * thick 
buckling_SF = bucklingcritmoment / M

max_bending_aly = M*R / (pi/4 * (R^4 - (R-thick)^4)) ;
max_shear_aly = 2*V / (pi * (R^2 - (R-thick)^2)^2);

sfbend = min(Xt, Xc) / max_bending_aly
sfshear = S / max_shear_aly

alumarea = pi*(R^2 - (R-thick)^2);
weight = alumarea * rho


function [sigmacr] = sigmacrit(E, v, r, t)
    sigmacr = E / sqrt(3*(1-v^2)) * t / r;
end

