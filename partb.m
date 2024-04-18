clear
close
%% Variables
M = 15e6;
V = 1.5e6;
R = 3;

Ex = 142e9; 
Ey = 11.2e9;
Gxy = 5;
vxy = 0.3;

Xt = 2200e6;
Xc = 1800e6;
Yt = 70e6;
Yc = 300e6;
S = 100e6;
t = 0.135e-3;
rho = 1610;


%% Skin Layups
%%Skin layups

% layup1 = zeros([1,180]) + 45;
% layup1(1:3) = [45, -45, 45];
% layup1(97:165) = 0;
% layup1(35:55) = 45;

l1 = [45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45];
l2 = [30, -30, 30, -30, 30, -30, 30, -30];
l3 = [60, -60, 60, -60, 60, -60, 60, -60];
l4 = [0, 90, 0, 90, 0, 90];
l5 = [l1, l1, l2, l3, l4];
l6 = flip(l5, 2);
layup1 = [l5, l6];

% layup1 = [45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 45, -45, 45, -45, 30, -30, 30, -30, 45, -45, 30, -30, 30, -30, 45, -45, 60, -60, 60, -60, 15, -15, 0, 90, 0, 90, 45, -45, 90, 0];
% layup1 = [45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45,45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45];

thicknesses1 = thicknessesgen(layup1, t);
t_skin = numel(thicknesses1) * t ;

thicknessm = numel(thicknesses1) * t * 1000;
skinarea = pi*(R^2 - (R-sum(thicknesses1))^2);


Q_lam1 = Qlam(Ex, Ey, vxy, Gxy);
Q_array1 = Qarray(layup1, Q_lam1);
ABD1 = ABD_matrix(layup1, thicknesses1, Q_array1);
A1 = ABD1(1:3, 1:3);
D1 = ABD1(4:6,4:6);


b = 0.8;
mof = 1.1;

I_z = pi/4 * (R^4 - (R - t_skin)^4);
stress_bottom = -M * R / I_z;

F_equiv = stress_bottom * t_skin

FI1 = safetyfact(2 * F_equiv, ABD1, thicknesses1, layup1, Q_lam1, Xt,...
    Xc, Yt, Yc, vxy, mof, Ex, S);

AR = linspace(0.3, 3, 100);

for m = 1:15
    for i = 1:numel(AR)
        ARspec = AR(i);
        a = b*ARspec;
        N0_top(m,i) = platebucklingccss(D1, ARspec, a, m);
    end
end
load = min(N0_top,[],1);
plot(AR, load)
hold on 
F = zeros(size(AR)) - F_equiv;
plot(AR, F)
title('AR vs Critical Buckling Load')

buckling_load = min(N0_top,[],"all");

a = 0.5;
b = 0.8;
AR = a/b;
m=1;
Ncrit = platebucklingccss(D1, AR, a, m) % N/m
FI1

bucklingSF =   abs(Ncrit / F_equiv)   % N/m / N/m
 
weight = skinarea * rho




%% Functions

function [sf] = safetyfact(F, ABD, thicknesses, layup, Q_lam, Xt, Xc, Yt, Yc, v12, mof, E1, S)
    load = [F;0;0;0;0;0];
    midplane_strain = linsolve(ABD, load);

    [strains_glob, strains_princ, stresses_glob, stresses_princ...
    ] = ply_strains(midplane_strain, thicknesses, layup, Q_lam);


    for k=1:size(stresses_princ,2)
            sigma1 = stresses_princ(1,k); sigma2 = stresses_princ(2,k);
            sigma3 = 0; sigma12 = stresses_princ(3,k);
            ff(k) = puck_ff(sigma1, sigma2, sigma3, Xt, Xc, v12, mof, E1);
            iff(k) = puck_iff(sigma2, sigma12, Yt, Yc, S);
    end

    sf = max(max(ff), max(iff));

end

function [sigmacr] = sigmacrit(E, v, r, t)
    sigmacr = E / sqrt(3*(1-v^2)) * t / r;
end

function [M0] = platebucklingmoment
    lambda = a/b * (D22 / D11)^0.25;
    K1 = 0.47 * pi^2 * b^2;
    K2 = m^2/lambda^2 + 2*(D12 + 2*D66)/(sqrt(D11*D22)) + lambda^2/m^2;
    K3 = m^2/lambda^2 + 2*(D12 + 2*D66)/(sqrt(D11*D22)) + 16*lambda^2/m^2;

    K = K1 * sqrt(K2 * K3);

    M0 = pi^2 * sqrt(D11*D22) * K / b^2;
    
end

function [N0] = platebucklingsscc(D, AR, a, m)
    D11 = D(1,1);
    D12 = D(1,2);
    D66 = D(3,3);
    D22 = D(2,2);

    b = a/AR;
    
    lambda = a/b * (D22 / D11)^0.25;
    K = m^2/lambda^2 + 2*(D12 + 2*D66)/(sqrt(D11*D22)) + 16/3*lambda^2/m^2;

    N0 = pi^2 / b^2 * sqrt(D11 * D22) * K;

end

function [N0] = platebucklingccss(D, AR, a, m)
    D11 = D(1,1);
    D12 = D(1,2);
    D66 = D(3,3);
    D22 = D(2,2);
    b = a/AR;

    lambda = a/b * (D22 / D11)^0.25;

    if lambda < 1.662
        K = m^2/lambda^2 + 2*(D12 + 2*D66)/(sqrt(D11*D22)) + 16/3*lambda^2/m^2;
    else
        k1 = (m^4 +8*m^2 + 1) / (lambda^2 * (m^2 + 1));
        k2 = 2*(D12 + 2*D66)/(sqrt(D11*D22));
        k3 = lambda^2 / (m^2 + 1);
        K = k1 + k2 + k3;
    end

    N0 = pi^2 / b^2 * sqrt(D11 * D22) * K;
end

function [N0] = platebucklingssuniax(D, AR, a, m)
    D11 = D(1,1);
    D12 = D(1,2);
    D66 = D(3,3);
    D22 = D(2,2);

    num = pi^2 * (D11*m^4 +2*(D12 + 2*D66)*m^2*AR^2 + D22*AR^4);
    denom = a^2*m^2;
    
    N0 = num / denom;
end



function [thick] = thicknessesgen(layup, t)
    thick = zeros(1,numel(layup));
    for i = 1:numel(layup)
        thick(i) = t;
    end
end

function [z] = zk(thickness)
    h = abs( sum(thickness)/2 );
    z = -h;
    for k = 1:numel(thickness)
        z = [z z(end)+thickness(k)];

    end
end

function [Q_lam] = Qlam(E1, E2, v12, G12)
    v21 = v12 * E2 / E1;

    
    Q = 1 - v12*v21;
    Q11 = E1 / Q;
    Q22 = E2 / Q;
    Q12 = v12 * E2 / Q;
    Q66 = G12;

    Q_lam = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];
end

function [Q] = Q_matrix(theta_angle, Q_lam)
    m = cos(deg2rad(theta_angle));
    n = sin(deg2rad(theta_angle));
    
    Q11 = Q_lam(1,1);
    Q12 = Q_lam(1,2);
    Q22 = Q_lam(2,2);
    Q66 = Q_lam(3,3);

    Qxx = Q11*m^4 + 2 * (Q12 + 2*Q66)*(m^2)*(n^2) + Q22*n^4;
    Qxy = (Q11 + Q22 - 4*Q66)*(m^2)*(n^2) + Q12*(m^4 + n^4);
    Qyy = Q11*n^4 + 2*(Q12 + 2*Q66)*(m^2)*(n^2) + Q22*m^4;
    Qxs = (Q11-Q12-2*Q66)*n*m^3 + (Q12-Q22+2*Q66)*n^3*m;
    Qys = (Q11-Q12-2*Q66)*m*n^3 + (Q12-Q22+2*Q66)*m^3*n;
    Qss = (Q11 + Q22 - 2*Q12 - 2*Q66)*(m^2)*(n^2) + Q66*(n^4 + m^4);
    
    Q = zeros(3);
    Q(1,1) = Qxx;
    Q(1,2) = Qxy;
    Q(2,1) = Qxy;
    Q(2,2) = Qyy;
    Q(3,3) = Qss;
    Q(3,1) = Qxs;
    Q(1,3) = Qxs;
    Q(3,2) = Qys;
    Q(2,3) = Qys;
    
end

function [Q_array] = Qarray(theta, Q_lam)
    % k = numel(theta);
    Q_array = zeros(3,3,8);
    
    for k = 1:numel(theta)
        angle = theta(k);
        Q_array(:,:,k) = Q_matrix(angle, Q_lam);
    end
end

function [ABD] = ABD_matrix(theta, thickness, Q_array)
    tol = 1e-5;
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    z = zk(thickness);

    for k = 1:numel(theta)
        
        Q = Q_array(:,:,k);
    
        A = A + Q .* thickness(k);
        B = B + (1/2) * Q .* ( z(k+1)^2 - z(k)^2 );
        D = D + (1/3) * Q .* ( (z(k+1))^3 - (z(k))^3 );

    end

    ABD = [A B; B D]; 

    ABD(ABD<tol)=0;
end

function [global_prop] = glob_prop(A, thickness_array)
    Ex = ( A(1,1)*A(2,2) - A(1,2)^2) / (sum(thickness_array)*A(2,2));
    Ey = ( A(1,1)*A(2,2) - A(1,2)^2) / (sum(thickness_array)*A(2,2));
    vxy = A(1,2) / A(2,2);
    vyx = A(1,2) / A(1,1);
    Gxy = A(3,3) * sum(thickness_array);
    global_prop = [Ex, Ey, vxy, vyx, Gxy];
end

function [prin_args] = direc_matrix(glob_args, theta)
    m = cos(deg2rad(theta));
    n = sin(deg2rad(theta));
    direc = [m^2, n^2, m*n; n^2, m^2, -m*n; -2*m*n, 2*m*n, (m^2 - n^2)];

    prin_args = direc * glob_args;
end

function [z_list, angle_at_z] = point_comp(theta, thickness)
    angle_at_z=[theta; theta];
    angle_at_z=angle_at_z(:)';
   
    z = zk(thickness);
    z_list=[z; z];
    z_list = z_list(:)';
    z_list = z_list(2:end-1);

end

function [strains_glob, strains_princ, stresses_glob, stresses_princ...
    ] = ply_strains(midplane_strain, thickness, theta, Q_lam)
    
    [z_list, angle_at_z] = point_comp(theta, thickness);

    %Strain Calc
    epsilonxx = midplane_strain(1) + z_list*midplane_strain(4);
    epsilonyy = midplane_strain(2) + z_list*midplane_strain(5);
    gammaxy = midplane_strain(3) + z_list*midplane_strain(6);

    strains_glob = [epsilonxx; epsilonyy; gammaxy];

    strains_princ = zeros(size(strains_glob));

    for i=1:numel(z_list)
        theta = angle_at_z(i);
        strains = strains_glob(:,i);
        strains_princ(:,i) = direc_matrix(strains, theta);
    end
    
    %Stress Calc
    stresses_glob = zeros(size(strains_glob));

    for i = 1:numel(z_list)     
        stresses_glob(:,i) = Q_matrix(angle_at_z(i), Q_lam) * strains_glob(:,i);
    end

    stresses_princ = zeros(size(stresses_glob));
    for k=1:numel(z_list)
        stresses_princ(:,k) = Q_lam * strains_princ(:,k);
    end
  
end
function [failure] = puck_ff(sigma1, sigma2, sigma3, Xt, Xc, v12, mof, E1)
    
    % Assume fibre poisson is 0.2, E_f = 225
    vf = 0.2;
    Ef = 225e9;

    %Initialise logical representing whether its failed. 
    failure = false;

    if sigma1 > 0
        R = Xt;
    else
        R = Xc;
    end

    a = (1/R);
    b = sigma1 - (v12 - vf*mof* E1/Ef)*(sigma2+sigma3);
    

    failure = a * b;

end

function [failure] = puck_iff(sigma2, sigma12, Yt, Yc, S)
    failure = false;
    p12p = 0.3;
    p12n = 0.2;

    %Mode A
    if sigma2 > 0
        p12p = 0.3;
        a = (sigma12/S)^2;
        b = (1-p12p*Yt/S)^2;
        c = (sigma2/Yt)^2;
        d = p12p * sigma2 / S;

        f = sqrt(a + b * c) + d;

    %Mode B and C 
    elseif sigma2 < 0
        sigma23a = S / (2 * p12n) * (sqrt(1+2*p12n*Yc/S)-1);
        p23n = p12n * sigma23a / S;
        sigma12c = S* sqrt(1 + 2*p23n);

        if abs(sigma2/sigma12) >= 0 && abs(sigma2/sigma12) <= sigma23a / abs(sigma12c)
            %Mode B
            f = 1/S * (sqrt(sigma12^2 + (p12n * sigma2)^2) + p12n*sigma2);
        else
            %Mode C
            f = ((sigma12 / (2*(1 + p23n*S)))^2 + ...
                (sigma2 / Yc)^2) * Yc / (-sigma2);
        end
    end

    failure = f;
end