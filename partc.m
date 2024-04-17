clear 

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

thick = 2.5e-3;

%% Skin stiffend panel buckling Analysis
a = 0.5;
b_spacing = 1.5;
AR = a/b_spacing;

l1 = [45, -45, 45, -45, 45, -45, 45, -45, 45, -45];
l2 = [30, -30, 30, -30, 30, -45];
l3 = [60, -60];
l4 = [0, 90, 0, 90, 0, 90, 0, 90];
l5 = [l1, l2, l3, l4];
l6 = flip(l5, 2);
skinlayup = [l5, l6];

t11 = [45, -45, 45, -45, 45, -45];
t12 = [60, -60,  30, -30];
t13 = [0, 90, 0, 90, 0, 90];
t14 = [t11, t12, t13];
t15 = flip(t14);
tsection1 = [t14, t15];

t21 = [45, -45, 45, -45];
t22 = [60, -60];
t23 = [0, 90, 0, 90, 0, 90];
t24 = [t21, t22, t23];
t25 = flip(t24);
tsection2 = [t24, t25];

section1_b = 80e-3;
section2_b = 110e-3;


%Thickness Generation
thicknessesskin= thicknessesgen(skinlayup, t);
t_skin = sum(thicknessesskin)
thicknessessec1= thicknessesgen(tsection1, t);
t_sec1= sum(thicknessessec1)
thicknessessec2= thicknessesgen(tsection2, t);
t_sec2 = sum(thicknessessec2)

% Total Force
sf = 2;
F_equiv = M * R * t_skin / (pi/4 * (R^4 - (R - t_skin)^4))
Ftot = F_equiv * sf;


%Skin
Q_lamskin = Qlam(Ex, Ey, vxy, Gxy);
Q_arrayskin = Qarray(skinlayup, Q_lamskin);
ABDskin = ABD_matrix(skinlayup, thicknessesskin, Q_arrayskin);
Askin = ABDskin(1:3,1:3);
Dskin = ABDskin(4:6,4:6);

m=4;
Pcr = platebucklingccss(Dskin, AR, a, m);

skin_bucklingSF = Pcr / F_equiv


bskin = a / 2*(1+2 * (a+Askin(2,1)/Askin(1,1)) * (1 - Pcr/Ftot) * ...
    (Askin(1,1)/(Askin(1,1)+3*Askin(2,2))))
EA_axialskin = EA(t_skin, Askin, bskin);

abdskin = ABDskin^-1;
EIskin = EI(t_skin, abdskin(4:6, 4:6));



%Section 1 : adjacent to skin
Q_lamsec1 = Qlam(Ex, Ey, vxy, Gxy);
Q_arraysec1 = Qarray(tsection1, Q_lamsec1);
ABDsec1 = ABD_matrix(tsection1, thicknessessec1, Q_arraysec1);
EA_sec1 = EA(t_sec1, ABDsec1(1:3,1:3), section1_b);
abdsec1 = ABDsec1^-1 ;
EIsec1 = EI(t_sec1, abdsec1(4:6, 4:6));

%Section 2 : vertical to skin
Q_lamsec2 = Qlam(Ex, Ey, vxy, Gxy);
Q_arraysec2 = Qarray(tsection2, Q_lamsec2);
ABDsec2 = ABD_matrix(tsection2, thicknessessec2, Q_arraysec2);
EA_sec2 = EA(t_sec2, ABDsec2(1:3,1:3), section2_b);
abdsec2 = ABDsec2^-1 ;
EIsec2 = EI(t_sec2, abdsec2(4:6, 4:6));

%Force Distributions 
Fskin = Ftot * EA_axialskin / (EA_sec1 + EA_sec2);
Fsec1 = Ftot * EA_sec1 / (EA_axialskin + EA_sec2);
Fsec2 = Ftot * EA_sec2 / (EA_axialskin + EA_sec1);

%Puck failure checking
mof = 1.1;

FIskin = safetyfact(-Fskin / (bskin*2), ABDskin, thicknessesskin, skinlayup, Q_lamskin, Xt, Xc, Yt, Yc, vxy, mof, Ex, S)
FIsec1 = safetyfact(-Fsec1 / section1_b, ABDsec1, thicknessessec1, tsection1, Q_lamsec1, Xt, Xc, Yt, Yc, vxy, mof, Ex, S)
FIsec2 = safetyfact(-Fsec2 / section2_b, ABDsec2, thicknessessec2, tsection2, Q_lamsec2, Xt, Xc, Yt, Yc, vxy, mof, Ex, S)


EI_equiv = EIskin + EIsec1 + EIsec2;

Pcrstiff = 7.56 * pi^2 * EI_equiv / a^2;

bucklingSF = Pcrstiff / Ftot


% Crippling Calc
D66crip = ABDsec1(6,6);
Nxcrip = 12 / D66crip / (0.5 * section1_b)^2;
sigmacrip = Nxcrip / (0.5 * section1_b);

ratio = 1.63 / (0.5 * section1_b / t_sec1)^0.717;

sigmaultc = sigmacrip / ratio; 
%Re run puck criteria with indexed compressive strength
FIcrip = safetyfact(-Fsec1 / section1_b, ABDsec1, thicknessessec1, tsection1, Q_lamsec1, Xt, sigmaultc, Yt, Yc, vxy, mof, Ex, S)


% Mass Calc

skinarea = pi*(R^2 - (R-t_skin)^2);
stiffarea = section1_b * t_sec1 + section2_b * t_sec2;
nstiff = 8;
totarea = skinarea + nstiff*stiffarea;
weight = totarea * rho



function [sf] = safetyfact(F, ABD, thicknesses, layup, Q_lam, Xt, Xc, Yt, Yc, v12, mof, E1, S)
    load = [F;0;0;0;0;0];
    midplane_strain = linsolve(ABD, load);

    [strains_glob, strains_princ, stresses_glob, stresses_princ...
    ] = ply_strains(midplane_strain, thicknesses, layup, Q_lam);

    for k=1:size(stresses_princ,2)
            sigma1 = stresses_princ(1,k); 
            sigma2 = stresses_princ(2,k);
            sigma3 = 0; 
            sigma12 = stresses_princ(3,k);
            ff(k) = puck_ff(sigma1, sigma2, sigma3, Xt, Xc, v12, mof, E1);
            iff(k) = puck_iff(sigma2, sigma12, Yt, Yc, S);
    end

    sf = max(max(ff), max(iff));

end

function [EAi] = EA(t, A, b)
    EAi = 1/t * (A(1,1) - A(2,1)^2 / A(2,2)) * b * t;
end

function [EIi] = EI(t, d)
    d11 = d(1,1);
    EIi = 12 / ( t^3 * d11) ;
end

% Buckling load of skin
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