% Clear workspace and close figures
clc; clear; close all;

% ===================== Environment Conditions =====================
rho = 0.0008754; % Air density (slug/ft^3)
T = 400;         % Temperature (°R)
mu = 3.025E-7;   % Dynamic viscosity (lb*s/ft^2)

% ===================== Geometry Parameters =====================
% Aircraft weight
W = 98000; % Weight (lb)

% Wing geometry
b_wing = 97;                % Wingspan (ft)
S_ref = 580.545 * 2;        % Reference area (ft^2)
t_c_wing = 0.1374;    % Thickness-to-chord ratio
swept_angle_wing = 6;       % Sweep angle (degrees)
sigma_wing = 0.26;          % Taper ratio
C_r_0_wing = 19;            % Root chord length (ft)

% Fuselage geometry
L_f = 92;                       % Fuselage length (ft)
D_f = 11;                       % Fuselage diameter (ft)
S_wet_f = 0.8 * pi * D_f * L_f; % Wetted area (ft^2)

% Horizontal tail geometry
S_ref_ht = 261;            % Horizontal tail area (ft^2)
t_c_ht = 0.09;             % Thickness-to-chord ratio
swept_angle_ht = 31.6;     % Sweep angle (degrees)
sigma_ht = 0.35;           % Taper ratio
C_r_ht = 11.1;             % Root chord length (ft)

% Vertical tail geometry
S_ref_vt = 161;            % Vertical tail area (ft^2)
t_c_vt = 0.09;             % Thickness-to-chord ratio
swept_angle_vt = 43.5;     % Sweep angle (degrees)
sigma_vt = 0.8;            % Taper ratio
C_r_vt = 15.5;             % Root chord length (ft)

% Pylon geometry
S_wet_p = 117;          % Wetted area (ft^2)
t_c_p = 0.06;           % Thickness-to-chord ratio
swept_angle_p = 0;      % Sweep angle (degrees)
sigma_p = 1;            % Taper ratio
c_p = 16.2;             % Chord length (ft)

% Nacelle geometry
S_wet_n = 455;          % Wetted area (ft^2)
Ln_Dn = 5;              % Length-to-diameter ratio
L_n = 16.8;             % Nacelle length (ft)

% ===================== Arrays for Graphing Results =====================
D_p = [];       % Parasite drag
D_i = [];       % Induced drag
D_total = [];   % Total drag
L_D = [];       % Lift-to-drag ratio

% ===================== Loop over Velocity Range =====================
for v = 230:880 % True airspeed (ft/s)
    % ===== (a) Parasite Drag Coefficient Calculation =====
    R_L = (rho * v) / mu;        % Reynolds number over characteristic length
    q = 0.5 * rho * v^2;         % Dynamic pressure (lb/ft^2)
    
    % ----- Wing Parasite Drag -----
    C_t_wing = sigma_wing * C_r_0_wing;      % Tip chord
    C_r_e_wing = 17.4056;                    % Root chord from model (ft)
    MAC_wing = (2/3) * (C_r_e_wing + C_t_wing - (C_r_e_wing * C_t_wing)/(C_r_e_wing + C_t_wing)); % Mean aerodynamic chord (ft)
    Re_wing = R_L * MAC_wing;                % Reynolds number (wing)
    C_f_wing = 0.0798 * Re_wing^-0.195;      % Skin friction coefficient (wing)
    Z = ((2 - 0.5^2) * cosd(swept_angle_wing)) / sqrt(1 - 0.5^2 * cosd(swept_angle_wing)^2); 
    K_wing = 1 + Z * t_c_wing + 100 * t_c_wing^4;  % Form factor
    S_wet_wing = 2*(498.95 + 496.3);         % Wetted area from model (ft^2)
    f_wing = K_wing * C_f_wing * S_wet_wing; % Equivalent profile drag area
    C_Dp_wing = f_wing / S_ref;              % Parasite drag coefficient (wing)

    % ----- Fuselage Parasite Drag -----
    Re_f = R_L * L_f;                        % Reynolds number (fuselage)
    C_f_f = 0.0798 * Re_f^(-0.195);          % Skin friction coefficient (fuselage)
    Lf_Df = L_f / D_f;                       % Length-to-diameter ratio
    K_f = 2.29 - 0.353 * Lf_Df + 0.038 * Lf_Df^2 - 1.48E-03 * Lf_Df^3; % Form factor (fuselage)
    f_f = K_f * C_f_f * S_wet_f;             % Equivalent profile drag area
    C_Dp_fuselage = f_f / S_ref;             % Parasite drag coefficient (fuselage)
    
    % ----- Horizontal Tail Parasite Drag -----
    MAC_ht = (2/3) * C_r_ht * (1 + sigma_ht - sigma_ht / (1 + sigma_ht)); % Mean aerodynamic chord (ft)
    Re_ht = R_L * MAC_ht;                   % Reynolds number (horizontal tail)
    C_f_ht = 0.0798 * Re_ht^-0.195;         % Skin friction coefficient (horizontal tail)
    Z_ht = ((2 - 0.5^2) * cosd(swept_angle_ht)) / sqrt(1 - 0.5^2 * cosd(swept_angle_ht)^2); 
    K_ht = 1 + Z_ht * t_c_ht + 100 * t_c_ht^4; % Form factor
    S_wet_ht = 2 * 1.02 * S_ref_ht;         % Wetted area (ft^2)
    f_ht = K_ht * C_f_ht * S_wet_ht;        % Equivalent profile drag area
    C_Dp_horizontal_tail = f_ht / S_ref;    % Parasite drag coefficient (horizontal tail)

    % ----- Vertical Tail Parasite Drag -----
    MAC_vt = (2/3) * C_r_vt * (1 + sigma_vt - sigma_vt / (1 + sigma_vt)); % Mean aerodynamic chord (ft)
    Re_vt = R_L * MAC_vt;                   % Reynolds number (vertical tail)
    C_f_vt = 0.0798 * Re_vt^-0.195;         % Skin friction coefficient (vertical tail)
    Z_vt = ((2 - 0.5^2) * cosd(swept_angle_vt)) / sqrt(1 - 0.5^2 * cosd(swept_angle_vt)^2); 
    K_vt = 1 + Z_vt * t_c_vt + 100 * t_c_vt^4; % Form factor
    S_wet_vt = 2 * 1.02 * S_ref_vt;         % Wetted area (ft^2)
    f_vt = K_vt * C_f_vt * S_wet_vt;        % Equivalent profile drag area
    C_Dp_vertical_tail = f_vt / S_ref;      % Parasite drag coefficient (vertical tail)

    % ----- Nacelle Parasite Drag -----
    Re_n = R_L * L_n;                     % Reynolds number (nacelle)
    C_f_n = 0.0798 * Re_n^-0.195;         % Skin friction coefficient (nacelle)
    K_n = 2.29 + -0.353*Ln_Dn + 0.038*Ln_Dn^2 + -1.48E-03*Ln_Dn^3;         % Form factor based on length-to-diameter ratio
    f_n = K_n * C_f_n * S_wet_n;          % Equivalent profile drag area for nacelle
    C_Dp_nacelles = f_n / S_ref;          % Parasite drag coefficient (nacelle)

    % ----- Pylon Parasite Drag -----

    Re_p = R_L * c_p;                     % Reynolds number (Pylon)
    C_f_p = 0.0798*Re_p^(-0.195);         % Skin friction coefficient (Pylon)
    Z_p = ( (2-0.5^2)*cos(deg2rad(swept_angle_p)) ) / sqrt(1 - 0.5^2*(cos(deg2rad(swept_angle_p)))^2);
    K_p = 1 + Z_p*t_c_p + 100*t_c_p^4;    % Form factor
    f_p = K_p * C_f_p * S_wet_p;          % Equivalent profile drag area for Pylon
    C_Dp_pylon = f_p/S_ref;               % Parasite drag coefficient (Pylon)

    % ===== Total Parasite Drag Coefficient =====
    C_Dp_total = 1.1 * (C_Dp_wing + C_Dp_fuselage + C_Dp_horizontal_tail + C_Dp_vertical_tail + C_Dp_nacelles + C_Dp_pylon);

    % ===== (b) Induced Drag Coefficient Calculation =====
    C_L = W / (q * S_ref);                 % Lift coefficient
    AR = b_wing^2 / S_ref;                 % Aspect ratio
    
    % ----- Oswald Efficiency Factor (e) -----
    if 0.0100 < C_Dp_total && C_Dp_total < 0.0150
        e1 = 0.969 - 0.0117 * AR + 1.85E-04 * AR^2;
        e15 = 0.975 - 0.0184 * AR + 3.7E-04 * AR^2;
        e = (e15 - e1) * (C_Dp_total - 0.0100) / (0.0150 - 0.0100) + e1;
        
    elseif 0.0150 < C_Dp_total && C_Dp_total < 0.0200
        e15 = 0.975 - 0.0184 * AR + 3.7E-04 * AR^2;
        e2 = 0.97 - 0.0226 * AR + 4.4E-04 * AR^2;
        e = (e2 - e15) * (C_Dp_total - 0.0150) / (0.0200 - 0.0150) + e15;
        
    elseif 0.0200 < C_Dp_total && C_Dp_total < 0.0250
        e2 = 0.97 - 0.0226 * AR + 4.4E-04 * AR^2;
        e25 = 0.958 - 0.0247 * AR + 4.07E-04 * AR^2;
        e = (e25 - e2) * (C_Dp_total - 0.0200) / (0.0250 - 0.0200) + e2;
        
    else
        disp('ERROR IN Oswald efficiency e');
    end


    % Induced drag coefficient
    C_Di = C_L^2 / (pi * AR * e);

    % ===== Forces =====
    D_p(end+1) = q * S_ref * C_Dp_total;       % Parasite drag
    D_i(end+1) = q * S_ref * C_Di;             % Induced drag
    C_D_total = C_Di + C_Dp_total;             % Total drag coefficient
    D_total(end+1) = q * S_ref * C_D_total;    % Total drag
    L_D(end+1) = C_L / C_D_total;              % Lift-to-drag ratio
end

% ===================== Plotting Results =====================
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
