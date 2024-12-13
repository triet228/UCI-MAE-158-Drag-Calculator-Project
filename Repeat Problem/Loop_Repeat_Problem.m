clc;clear;close all;
% Given Parameter
% Environment condition
rho = 0.0008754; % slug/ft^3
T = 400; % °R
mu = 3.025E-7; % lb∗s / ft^2

% Geometry parameter
W = 98000;
% Wing ----------------------------Change later
b_wing = 93.2; % ft
S_ref = 1000; % ft^2
t_c_wing = 0.106;
swept_angle_wing = 24.5;
sigma_wing = 0.2;
C_r_0_wing = 17.8; % ft
expose = 1 - 0.17;

% Fuselage ----------------------------Change later
L_f = 107;
D_f = 11.5;
S_wet_f = 3280;

% Horizontal tail
S_ref_ht = 261;
t_c_ht = 0.09;
swept_angle_ht = 31.6;
sigma_ht = 0.35;
C_r_ht = 11.1;

% Vertical tail
S_ref_vt = 161;
t_c_vt = 0.09;
swept_angle_vt = 43.5;
sigma_vt = 0.8;
C_r_vt = 15.5;

% Pylons
S_wet_p = 117;
t_c_p = 0.06;
swept_angle_p = 0;
sigma_p = 1;
c_p = 16.2;

% Nacelles
S_wet_n = 455;
Ln_Dn = 5;
L_n = 16.8;


% Prepare list for Graph
D_p = [];
D_i = [];
D_total = [];
L_D = [];


for v = 230:880 % True air speed

% $$$$$$$$$$$$$$$$ (a) Parasite drag coefficient $$$$$$$$$$$$$$$$$$$$$$$$

% Reynolds number over characteristic length (Re/L)
R_L = (rho*v) /mu;

% Dynamic pressure
q = 0.5*rho*v^2;

% $$$$$$$$$$$$$$$$ Wing $$$$$$$$$$$$$$$$

% Tip cord
C_t_wing = sigma_wing*C_r_0_wing;
% Root cord from Solidworks file (Wing_MAC.SLDPRT)----------------------------Change later
C_r_e_wing = 16.04291846; % ft
% Mean aerodynamic chord   ----------------------------Change later
MAC_wing = (2/3) * (C_r_e_wing + C_t_wing - (C_r_e_wing * C_t_wing)/(C_r_e_wing + C_t_wing) ); % ft

% Reynolds number
Re_wing = R_L * MAC_wing;

% Skin friction coefficient Cf
C_f_wing = 0.0798*Re_wing^-0.195;

% Form factor K
Z = ( (2-0.5^2)*cos(deg2rad(swept_angle_wing)) ) / sqrt(1 - 0.5^2*(cos(deg2rad(swept_angle_wing)))^2);
K_wing = 1 + Z*t_c_wing + 100*t_c_wing^4;

% Wetted area -----------------------------------------------------------Change later
S_wet_wing = 2*1.02*S_ref*expose;

% Equivalent profile drag area of the wing
f_wing = K_wing*C_f_wing*S_wet_wing;

% Parasite drag coefficient
C_Dp_wing = f_wing/S_ref;

% $$$$$$$$$$$$$$$$ Fuselage $$$$$$$$$$$$$$$$

% Reynolds number
Re_f = R_L * L_f;

% Skin friction coefficient Cf
C_f_f = 0.0798*Re_f^(-0.195);

% Form factor K
Lf_Df = L_f/D_f;
K_f = 2.29 + -0.353*Lf_Df + 0.038*Lf_Df^2 + -1.48E-03*Lf_Df^3;

% Equivalent profile drag area of the wing
f_f = K_f * C_f_f * S_wet_f;

% Parasite drag coefficient
C_Dp_fuselage = f_f/S_ref;


% $$$$$$$$$$$$$$$$ Horizontal Tail $$$$$$$$$$$$$$$$

% Mean aerodynamic chord   
MAC_ht = (2/3) * C_r_ht * (1 + sigma_ht - sigma_ht / ( 1 + sigma_ht)); % ft

% Reynolds number
Re_ht = R_L * MAC_ht;

% Skin friction coefficient Cf
C_f_ht = 0.0798*Re_ht^-0.195;

% Form factor K
Z_ht = ( (2-0.5^2)*cos(deg2rad(swept_angle_ht)) ) / sqrt(1 - 0.5^2*(cos(deg2rad(swept_angle_ht)))^2);
K_ht = 1 + Z_ht*t_c_ht + 100*t_c_ht^4;

% Wetted area 
S_wet_ht = 2*1.02*S_ref_ht;

% Equivalent profile drag area of the wing
f_ht = K_ht*C_f_ht*S_wet_ht;

% Parasite drag coefficient
C_Dp_horizontal_tail = f_ht/S_ref;

% $$$$$$$$$$$$$$$$ Vertical Tail $$$$$$$$$$$$$$$$

% Mean aerodynamic chord   
MAC_vt = (2/3) * C_r_vt * (1 + sigma_vt - sigma_vt / ( 1 + sigma_vt)); % ft

% Reynolds number
Re_vt = R_L * MAC_vt;

% Skin friction coefficient Cf
C_f_vt = 0.0798*Re_vt^-0.195;

% Form factor K
Z_vt = ( (2-0.5^2)*cos(deg2rad(swept_angle_vt)) ) / sqrt(1 - 0.5^2*(cos(deg2rad(swept_angle_vt)))^2);
K_vt = 1 + Z_vt*t_c_vt + 100*t_c_vt^4;

% Wetted area 
S_wet_vt = 2*1.02*S_ref_vt;

% Equivalent profile drag area of the wing
f_vt = K_vt*C_f_vt*S_wet_vt;

% Parasite drag coefficient
C_Dp_vertical_tail = f_vt/S_ref;


% $$$$$$$$$$$$$$$$ Nacelles $$$$$$$$$$$$$$$$

% Reynolds number
Re_n = R_L * L_n;

% Skin friction coefficient Cf
C_f_n = 0.0798*Re_n^(-0.195);

% Form factor K
K_n = 2.29 + -0.353*Ln_Dn + 0.038*Ln_Dn^2 + -1.48E-03*Ln_Dn^3;

% Equivalent profile drag area of the wing
f_n = K_n * C_f_n * S_wet_n;

% Parasite drag coefficient
C_Dp_nacelles = f_n/S_ref;


% $$$$$$$$$$$$$$$$ Pylon $$$$$$$$$$$$$$$$

% Reynolds number
Re_p = R_L * c_p;

% Skin friction coefficient Cf
C_f_p = 0.0798*Re_p^(-0.195);

% Form factor K
Z_p = ( (2-0.5^2)*cos(deg2rad(swept_angle_p)) ) / sqrt(1 - 0.5^2*(cos(deg2rad(swept_angle_p)))^2);
K_p = 1 + Z_p*t_c_p + 100*t_c_p^4;

% Equivalent profile drag area of the wing
f_p = K_p * C_f_p * S_wet_p;

% Parasite drag coefficient
C_Dp_pylon = f_p/S_ref;


% TOTAL Parasite drag coefficient
C_Dp_total = 1.1*(C_Dp_wing + C_Dp_fuselage + C_Dp_horizontal_tail + C_Dp_vertical_tail + C_Dp_nacelles + C_Dp_pylon);


% $$$$$$$$$$$$$$$$ (b) Induced drag coefficient $$$$$$$$$$$$$$$$$$$$$$$$

% Lift coefficient
C_L = W / (q*S_ref);

% Aspect Ratio
AR = b_wing^2 / S_ref;

% Oswald efficiency e
if 0.0100 < C_Dp_total & C_Dp_total < 0.0150
    e1 = 0.969 + -0.0117*AR + 1.85E-04*AR^2;
    e15 = 0.975 + -0.0184*AR + 3.7E-04*AR^2;
    e = (e15 - e1) * (C_Dp_total - 0.0100) / (0.0150  - 0.0100) + e1;
elseif 0.0150 < C_Dp_total & C_Dp_total < 0.0200
    e15 = 0.975 + -0.0184*AR + 3.7E-04*AR^2;
    e2 = 0.97 + -0.0226*AR + 4.4E-04*AR^2;
    e = (e2 - e15) * (C_Dp_total - 0.0150) / (0.0200  - 0.0150) + e15;
elseif 0.0200 < C_Dp_total & C_Dp_total < 0.0250
    e2 = 0.97 + -0.0226*AR + 4.4E-04*AR^2;
    e25 = 0.958 + -0.0247*AR + 4.07E-04*AR^2;
    e = (e25 - e2) * (C_Dp_total - 0.0200) / (0.0250  - 0.0200) + e2;
else
    disp('ERROR IN Oswald efficiency e');
end

% Induced drag coefficient
C_Di = C_L^2 / (pi*AR*e);

% $$$$$$$$$$$$$$$$ Forces $$$$$$$$$$$$$$$$$$$$$$$$
% Parasite drag
D_p(end+1) = q*S_ref*C_Dp_total;

% Induced drag
D_i(end+1) = q*S_ref*C_Di;

% Total Drag
C_D_total = C_Di + C_Dp_total;
D_total(end+1) = q*S_ref*C_D_total;

% Lift to Drag ratio
L_D(end+1) = C_L/C_D_total;

end

v_list = 230:880; % Velocity range for plotting

% Plot Parasite, Induced, and Total Drag vs Velocity
figure;
plot(v_list, D_p, v_list, D_i, v_list, D_total);
xlabel('Velocity (ft/s)');
ylabel('Drag (lb)');
legend('Parasite Drag', 'Induced Drag', 'Total Drag');
title('Drag vs. Velocity');

% Plot Lift-to-Drag Ratio vs Velocity
figure;
plot(v_list, L_D);
xlabel('Velocity (ft/s)');
ylabel('Lift-to-Drag Ratio');
title('Lift-to-Drag Ratio vs. Velocity');