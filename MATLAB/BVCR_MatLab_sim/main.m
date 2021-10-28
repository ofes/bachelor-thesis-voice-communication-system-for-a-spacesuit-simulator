%% MAIN: Self-sufficient energy distribution system performance estimation
% This simulation provides a performance estimation of the self-sufficient energy distribution system, which is part of the voice communication system of the OeWF Serenity spacesuit simulator. It can be simulated if the repeater radio infrastructure can operate at a certain mission location on Earth. All sections labeled [INPUT] can be manipulated. For further information please refer to the MATLAB documentation, look into the bachelor thesis about this simulation or contact the author.

% Organization:     OeWF (Austrian Space Forum)
% Author:           Omar Filip El Sendiouny
% Project:          Serenity BU-COMMs
% Date:             15.12.2020
% Version:          1 

clear all;
close all;
clc;

show_command_window_output = 1; % allow command window output
N_dp = 4; % number of decimal points for solar angle calculations (accuracy)

dfile ='command_window_output.txt'; % command window output can be found in the local folder
if(exist(dfile, 'file'))
    delete(dfile); 
end
diary(dfile)
diary on

%% [INPUT] Mission information:
% It it assumed that the missions start at the same time every day for the entire duration of the mission. The time of the mission start must be smaller than the time of mission end.

lat = 48.210; % latitude in (deg)
lon = 16.363; % longitude in (deg)
irradiation_GH = 3.25; % daily total of the global horizontal irradiation in (kWh / m^2) 
irradiation_DN = 2.9; % daily total of the direct normal irradiation in (kWh / m^2)
date(1) = datetime(2021, 4, 1); % mission start date in (Y, M, D)
date(2) = datetime(2021, 4, 30); % mission end date in (Y, M, D)
t_UTC(1) = {'5:00'}; % time of the daily mission start in (h)
t_UTC(2) = {'15:00'}; % time of the daily mission end in (h)
ALB = 0.25; % albedo in (1)
theta_A = 11.4; % average ambient temperature in (degrees C)

%% [INPUT] Photovoltaic generator:
% MPP ... maximim power point
% STC ... standard test conditions
% SC ... short-circuit
% OC ... open-circuit

N_C = 36; % number of PV cells in (1)
A_PV = 0.88; % PV generator area in (m^2)
I_MPP_STC = 9.11; % current at MPP for STC in (A)
U_MPP_STC = 21.41; % voltage at MPP for STC in (V)
I_SC_STC = 9.79; % SC current for STC in (A)
U_OC_STC = 24.27; % OC voltage for STC in (V)
TC_I_SC = 0.05; % temperature coefficient for I_SC in (% / degC)
TC_U_OC = -0.29; % temperature coefficient for U_OC in (% / degC)
NOCT = 45; % nominal operating cell temperature in (degC)
m = 1.19045; % ideality factor in (1)
N_PV = 2; % number of PV generators (if greater than 1, they are connected in parallel)

%% [IMPUT] Load (repeater radio system):
% After t_UTC(2) and before t_UTC(1) the repeater radio system is in standby. The sum of the duty cycles must be 100%.

U_L_min = 11.0; % minimum load voltage in (V)
U_L_max = 15.5; % maximum load voltage in (V)
I_T = 3; % current consumption when the load transmits in (A)
I_R = 0.7; % current consumption when the load receives in (A)
I_Stby = 0.7; % current consumption when the load is in standby in (A)
a_T = 20; % duty cycle when transmitting in (%) 
a_R = 20; % duty cycle when receiving in (%)
a_Stby = 60; % duty cycle when in standby in (%)
I_add = 0; % additional load current in (A)

%% [INPUT] LiFePo_4 battery
% The initial state of charge of the battery is usually the SOC for storage, which is between 50% and 60%.

SOC_init = 0.6; % initial state of charge in (1)
Q_nom = 50; % nominal charge in (Ah)
C_D = 1/3; % discharge rate for Q_nom in (h^(-1))

% Obtained from the discharging an charging experiments:
SOC = [1.00 0.95 0.90 0.85 0.80 0.75 0.70 0.65 0.60 0.55 0.50 0.45 0.40 0.35 0.30 0.25 0.20 0.15 0.10 0.05 0.00]; % state of charge in (1) of the battery
U_0_D = [13.877 13.315 13.319 13.322 13.320 13.306 13.248 13.196 13.181 13.117 13.174 13.172 13.169 13.111 13.087 13.012 12.931 12.837 12.814 12.350 11.371]; % open-circuit voltages when discharging the battery in (V)
U_0_C = [14.107 13.387 13.393 13.401 13.415 13.388 13.387 13.275 13.272 13.253 13.250 13.248 13.250 13.244 13.225 13.179 13.089 12.958 12.923 12.502 11.420]; % open-circuit voltages when charging the battery in (V)
R_e_D = [0.0120840740086789 0.0117154203104446 0.0117700588020157 0.0117132299485587 0.0116221338950336 0.0116467701188540 0.0116605558931997 0.0115638462618863 0.0116222251883952 0.0116974390402753 0.0116987272047564 0.0117465274490441 0.0117588908671468 0.0117566064015986 0.0117756658545291 0.0117553658832641 0.0117583775845481 0.0117768748688630 0.0117155094645921 0.0118414656950463 0.0120957986167104]; % electrolyte resistances when discharging in (Ohm)
R_e_C = [0.0121012406947893 0.0113356302935661 0.0114597641410590 0.0114014839745851 0.0114900762492556 0.0116012192457450 0.0112273867826839 0.0112396239841211 0.0111117629249463 0.0112231518757464 0.0113202903133204 0.0112805170265242 0.0114125559543423 0.0114503161497776 0.0114125257730806 0.0115293500897879 0.0115664251765242 0.0115250208156497 0.0117012683150631 0.0118226401014494 0.0120662511263975]; % electrolyte resistances when charging in (Ohm)

%% [INPUT] Solar charging controller
% The SCC only connects the PV generator when its voltage at the PV generator input terminal exceeds U_Bat + 5V. As soon as this condition is met the PV generator remains connected as long as the voltage at the PV generator input terminal is greater than U_Bat + 1V.

I_SCC = 0.020; % current consumption of the SCC in (A)
I_B_max = 10; % maximum battery current in (A)
P_MPP_max = 145; % maximum power from PV generator in (W)
eta_SCC = 0.98; % efficiency in (1)

%% [INPUT] Cables:
% Only one wire of the cable is to be described. For more information, please refer to the subsection 3.1.2 (Repeater radio infrastructure) of the bachelor thesis about this simulation.

% Input sytle: cables = [cable_A; cable_B; cable_C]
% Input style: cables(n,:) = [length cross_section_area specific_resistance temperature_coefficient]

% length in (m)
% cross_section_area in (mm^2)
% specific_resistance for 20 degC in (Ohm * mm^2 / m)
% temperature coefficientc for 20 degC in (degC^(-1))

cables = [0.9 4.0 0.01673 4.3 * 10^(-3); 0.5 4.0 0.01673 4.3 * 10^(-3); ...
    11.0 2.5 0.01673 4.3 * 10^(-3)];  

%% [OUTPUT] Radiation flux received by an inclined PV generator:
% To ensure better readability, for-loops were used instead of direct assignments for most vectors/matrices. This further helps with the maintenance of the program. However, it decreases the performance of the simulation (fix in version 1.1)

lat = round(lat, N_dp);
lon = round(lon, N_dp);

date_array = (date(1):date(2))'; % contains all mission dates
N_md = split(between(date(1), date(2), 'days'), 'days') + 1; % number of mission days

delta_t_S = 0.01; % resolution for solar time
N_dp_delta_t_S = 0; % number of decimal places of delta_t_S

while (floor(delta_t_S * 10^N_dp_delta_t_S) ~= delta_t_S * 10^N_dp_delta_t_S) % find the number of decimal places of delta_t_S
N_dp_delta_t_S = N_dp_delta_t_S + 1;
end

t_S = round(0:delta_t_S:24, N_dp_delta_t_S); % solar time
h_S = round((t_S - 12) * 15, N_dp); % solar hour angle

N_d = zeros(N_md, 1); % number of days since Jan. 1st
delta = zeros(N_md, 1); % Sun declination
beta = zeros(N_md, 1); % optimal inclination angle for the PV generator
h_rs = zeros(N_md, 2); % solar sunrise and sunset hour angles

for n = 1:N_md % calculating N_d, delta, and beta
    N_d(n) = 30.3 * (month(date_array(n)) - 1) + day(date_array(n));
    delta(n) = round(23.45 * sind((360 * (284 + N_d(n)) / 365)), N_dp);
    beta(n) = round(abs(lat - delta(n)), N_dp);
end

beta_mean = mean(beta); % mean value of beta will be used from now on

for n = 1:N_md % calculating h_rs
    if(abs(lat) < 90 - abs(delta(n)))
        h_rs(n,1) = round(- acosd(- tand(delta(n)) * tand(lat)), N_dp); % solar sunrise hour angle
        h_rs(n,2) = round(- h_rs(n,1), N_dp); % solar sunset hour angle
    else
        error('Solar sunrise and sunset angles cannot not be calculated!');
    end
end

t_rs = 12 + h_rs / 15; % solar sunrise and sunrise times

idx = h_S <= 0; % get index for solar hour angles from t_S = 0h to t_S = 12h
h_S_fhd = h_S(idx); % create associated vector for solar hour angles (first half of the day)

gamma_S = zeros(N_md, length(h_S_fhd)); % altitude of the Sun
alpha_S = zeros(N_md, length(h_S_fhd)); % azimuth of the Sun
theta = zeros(N_md, length(h_S_fhd)); % angle theta (solar incidence angle)

for n = 1:N_md % calculating gamma_S_deg from sunrise to t_S = 12h
    for o = 1:length(h_S_fhd)
        if(h_S_fhd(o) >= h_rs(n,1)) % if the Sun is above the horizon
            gamma_S(n,o) = round(asind(sind(lat) * sind(delta(n)) + cosd(lat) * cosd(delta(n)) * cosd(h_S_fhd(o))), N_dp); 
         else
            gamma_S(n,o) = NaN; % invalid angle (Sun is below the horizon)
        end
    end
end

for n = 1:N_md % calculating alpha_S_deg from sunrise to t_S = 12h
    for o = 1:length(h_S_fhd)
        if(h_S_fhd(o) >= h_rs(n,1) && h_S_fhd(o) ~= 0) % if the Sun is above the horizon before t_S = 12h
            alpha_S(n,o) = round(-acosd((sind(lat) * cosd(delta(n)) * cosd(h_S_fhd(o)) - cosd(lat) * sind(delta(n))) / cosd(gamma_S(n,o))), N_dp); 
        elseif(h_S_fhd(o) == 0 && gamma_S(n,o) ~= 90) % if t_S = 12h and the Sun is not at its zenith
            if(lat > delta(n))
                alpha_S(n,o) = round(0, N_dp);
            elseif(lat < delta(n))
                alpha_S(n,o) = round(180, N_dp);
            end
        else % the Sun is below the horizon or at its zenith
            alpha_S(n,o) = NaN; % invalid angle
        end  
    end
end

for n = 1:N_md % calculating theta
    for o = 1:length(h_S_fhd)
        if(h_S_fhd(o) >= h_rs(n,1)) % if the Sun is above the horizon
            if(gamma_S(n,t_S == 12) ~= 90) % if the Sun is not at its zenith
                tmp = sind(delta(n)) * sind(lat) * cosd(beta_mean) ...
                    - sind(delta(n)) * cosd(lat) * sind(beta_mean) * cosd(alpha_S(n, t_S == 12)) ...
                    + cosd(delta(n)) * cosd(lat) * cosd(beta_mean) * cosd(h_S_fhd(o)) ...
                    + cosd(delta(n)) * sind(lat) * sind(beta_mean) * cosd(alpha_S(n,o)) * cosd(h_S_fhd(o)) ...
                    + cosd(delta(n)) * sind(beta_mean) * sind(alpha_S(n, t_S == 12)) * cosd(h_S_fhd(o));
                if(tmp>1) % necessary because of rounding errors caused by matlab
                    tmp = 1;
                end 
                    theta(n,o) = round((acosd(tmp)), N_dp);      
            elseif(gamma_S(n,t_S == 12) == 90) % if the Sun is at its zenith
                tmp = sind(delta(n)) * sind(lat) ... 
                    + cosd(delta(n)) * cosd(lat) * cosd(h_S_fhd(o));
                theta(n,o) = round(acosd(tmp), N_dp);
            end
        else
            theta(n,o) = NaN; % invalid angle
        end
    end
end

for n = 1:N_md % because of the projection of the normal to A_PV onto plane earth?
    [~, idx(n)] = max(theta(n,:));
    for o = 1:length(theta)
        if(o < idx(n) || theta(n,o) >= 90.0) % if the Sun's rays do not hit A_PV
            theta(n,o) = NaN;
        end
    end
end

tmp_mat = gamma_S;
tmp_mat = fliplr(tmp_mat); % flip the temporary matrix from left to right
tmp_mat(:,1) = []; % delete first entry so it does not appear twice in the gamma_S matrix
gamma_S = [gamma_S tmp_mat]; % allowed because the Sun is symmetrical around t_S = 12h

tmp_mat = alpha_S;
tmp_mat = fliplr(tmp_mat); % flip the temporary matrix from left to right
tmp_mat(:,1) = []; % delete first entry so it does not appear twice in the alpha_S matrix
alpha_S = [alpha_S -tmp_mat]; % allowed because the Sun is symmetrical around t_S = 12h

tmp_mat = theta;
tmp_mat = fliplr(tmp_mat); % flip the temporary matrix from left to right
tmp_mat(:,1) = []; % delete first entry so it does not appear twice in the theta matrix
theta = [theta tmp_mat]; % allowed because the Sun is symmetrical around t_S = 12h

E_GHI = irradiation_GH * 10^(3) / 24; % global horizontal irradiance 
E_DNI = irradiation_DN * 10^(3) / 24; % direct normal irradiance

E_G = zeros(N_md, length(t_S)); % total irradiance received by an inclined PV generator
W_G = zeros(N_md, 1); % solar energy yield
tau_S = zeros(N_md, 1); % auxiliary variable
Phi_max = zeros(N_md, 1); % maximal daily radiation flux
Phi_G = zeros(N_md, length(t_S)); % radiation flux received by an inclined PV generator 

for n = 1:N_md % calculating E_G
    for o = 1:length(t_S)
        if(isnan(theta(n,o)) && ~isnan(gamma_S(n,o)))
            E_G(n,o) = (E_GHI - E_DNI * sind(gamma_S(n,o))) * (1 + cosd(beta_mean)) / 2; % PV generator receives DIFG
        elseif(~isnan(theta(n,o)) && ~isnan(gamma_S(n,o))) % PV generator receives DGI, DIFG and RGI
            E_G(n,o) = E_DNI * cosd(theta(n,o)) ...
                + (E_GHI - E_DNI * sind(gamma_S(n,o))) * (1 + cosd(beta_mean)) / 2 ...
                + E_GHI * ALB * (1 - cosd(beta_mean)) / 2;
        else
            E_G(n,o) = 0;
        end
    end
end

for n = 1:N_md % calculating W_G
    for o = 1:length(t_S)
        W_G(n) = W_G(n) + A_PV * E_G(n,o) * delta_t_S;
    end
end

for n = 1:N_md % calculating Phi_G
    for o = 1:length(t_S)
        tau_S(n) = (12 - t_rs(n,1)) / 3;
        Phi_max(n) = 0.997 / (tau_S(n) * sqrt(2 * pi)) * W_G(n);
        Phi_G(n,o) = Phi_max(n) * exp(- (t_S(o) - 12)^2/(2 * tau_S(n)^2));
    end
end
 
%% [OUTPUT] PV generator power output

e = 1.602176634 * 10^(-19); % elementary charge in (As)
k_B = 1.380649 * 10^(-23); % Bolzmann constant in (Ws/K)

theta_C = zeros(N_md, length(t_S)); % PV cell temperature
U_T = zeros(N_md, length(t_S)); % thermal voltage
I_Ph = zeros(N_md, length(t_S)); % photocurrent
I_Ph_irr_STC = zeros(N_md, length(t_S)); % photocurrent with E_STC
U_OC_0 = zeros(N_md, length(t_S)); % open-circuit voltage
I_S_0 = zeros(N_md, length(t_S)); % reverse saturation current

for n = 1:N_md % calculating the parameters initialized above
    for o = 1:length(t_S)
        theta_C(n,o) = theta_A + (NOCT - 20) * Phi_G(n,o) * A_PV / 800;
        U_T(n,o) = k_B * (theta_C(n,o) + 273.15)/e;
        I_Ph(n,o) = I_SC_STC * Phi_G(n,o) * A_PV * (1 + TC_I_SC / 100 * (theta_C(n,o) - 25)) / 1000;
        I_Ph_irr_STC(n,o) = I_SC_STC * (1 + TC_I_SC / 100 * (theta_C(n,o) - 25));
        U_OC_0(n,o) = U_OC_STC * (1 + TC_U_OC / 100 * (theta_C(n,o) - 25)) ...
            + m * N_C * U_T(n,o) * log(I_Ph(n,o) / I_Ph_irr_STC(n,o));
        
        if(U_OC_0(n,o) < 0) % so that I_S << I_Ph is fulfilled
            U_OC_0(n,o) = NaN;
            I_S_0(n,o) = NaN;
        else
            I_S_0(n,o) = I_Ph(n,o) * exp(- (U_OC_0(n,o)) / (m * N_C * U_T(n,o)));
        end
    end
end

[idxn, idxo] = find(U_OC_0 == max(max(U_OC_0))); % finding maximum value in U_OC_0 matrix (improve in version 1.1)
delta_U_PV = 0.01; % resolution for PV generator voltage
U_PV = 0:delta_U_PV:U_OC_0(idxn,idxo) + delta_U_PV; % PV generator voltage
I_PV = zeros(N_md, length(t_S), length(U_PV)); % PV generator current
P_PV = zeros(N_md, length(t_S), length(U_PV)); % PV generator power
P_MPP = zeros(N_md, length(t_S)); % PV generator power for MPP
U_MPP = zeros(N_md, length(t_S)); % PV generator voltage for MPP
I_MPP = zeros(N_md, length(t_S)); % PV generator current for MPP
W_G_el = zeros(N_md, 1); % electrical energy yield of the PV generator

for n = 1:N_md % caluculating PV generator current and MPP values as well as its daily electrical energy yield (3D matrix ... noice UwU)
    for o = 1:length(t_S)
        for p = 1:length(U_PV)
            I_PV(n,o,p) = N_PV * (I_Ph(n,o) - I_S_0(n,o) * (exp(U_PV(p)/(m * N_C * U_T(n,o))) - 1));
            
            if(I_PV(n,o,p) < 0)
                I_PV(n,o,p) = NaN;
            end
            
            P_PV(n,o,p) = I_PV(n,o,p) * U_PV(p);
            
            if(isnan(P_PV(n,o,p)))
                P_PV(n,o,p) = 0;
            end
        end
        
        [P_MPP(n,o), idx] = max(P_PV(n,o,:));
        U_MPP(n,o) = U_PV(idx);
        
        if(isnan(I_PV(n,o,idx)))
            I_MPP(n,o) = 0;
        else
            I_MPP(n,o) = I_PV(n,o,idx);
        end
    end
end

for n = 1:N_md % calculating W_G_el
    for o = 1:length(t_S)
        W_G_el(n) = W_G_el(n) + P_MPP(n,o) * delta_t_S; % electrical energy yield
    end
end

%% [OUTPUT] Load (repeater radio system)
% The load current is an average current, because an on-off step function could not be implemented.

dtv = datevec(datetime(t_UTC{1},'InputFormat','HH:mm')); % converitng time of daily mission start to floating point number
da = duration(dtv(:,4:end));
t_UTC_float(1) = hours(da);

dtv = datevec(datetime(t_UTC{2},'InputFormat','HH:mm')); % converitng time of daily mission end to floating point number
da = duration(dtv(:,4:end));
t_UTC_float(2) = hours(da);

Z_h = zeros(N_md, 1); % equation of time 
t_S_mission_se = zeros(N_md, 2); % solar mission start and end times
I_L = zeros(N_md, length(t_S)); % load current

for n = 1:N_md % calculating solar mission time
    Z_h(n) = 0.123 * cosd(360 * (88 + N_d(n)) / 365) - 0.167 * sind(720 * (10 + N_d(n)) / 365);
    t_S_mission_se(n, 1) = round(t_UTC_float(1) + Z_h(n) + lon / 15, N_dp_delta_t_S);
    t_S_mission_se(n, 2) = round(t_UTC_float(2) + Z_h(n) + lon / 15, N_dp_delta_t_S);
end

idx_mission_se = zeros(N_md, 2); % indices for daily mission start and end in t_S

for n = 1:N_md % find index of daily mission start and end in t_S
    [~,idx_mission_se(n,1)] = find(t_S == t_S_mission_se(n,1));
    [~,idx_mission_se(n,2)] = find(t_S == t_S_mission_se(n,2)); 
end

for n = 1:N_md % calculating I_L
    if(n == 1)
        for o = 1:length(t_S) % first day begins at mission start
            if(o < idx_mission_se(n,1))
                I_L(n,o) = 0;
            elseif(o >= idx_mission_se(n,1) && o <= idx_mission_se(n,2))
                I_L(n,o) = I_T * a_T / 100 + I_R * a_R / 100 + I_Stby * a_Stby / 100 + I_add;
            else
                I_L(n,o) = I_Stby + I_add;
            end
        end
    else
        for o = 1:length(t_S) % all other days
            if(o >= idx_mission_se(n,1) && o <= idx_mission_se(n,2))
                I_L(n,o) = I_T * a_T / 100 + I_R * a_R / 100 + I_Stby * a_Stby / 100 + I_add;
            else
                I_L(n,o) = I_Stby + I_add;
            end
        end
    end    
end

%% [OUTPUT] Cables losses

sz = size(cables);
R_wires = zeros(sz(1),1); % wire resistances

for p = 1:sz(1) % calculate all wire resistances
    R_wires(p) = cables(p,1) * cables(p,3) * (1 + cables(p,4) * (theta_A - 20)) ... 
        / cables(p,2);
end

%% [OUTPUT] Solar charging controller and LiFePO_4 battery

k_P = 1.05; % Peukert constant in (1)
I_nom = C_D * Q_nom; % nominal battery current

U_in = zeros(N_md, length(t_S)); % voltage at PV generator input terminal
U_L = zeros(N_md, length(t_S)); % voltage at repeater
Q_B = zeros(N_md, length(t_S)); % battery charge
I_B = zeros(N_md, length(t_S)); % battery current
U_B = zeros(N_md, length(t_S)); % battery voltage
SOC_curr = zeros(N_md, length(t_S)); % current state of charge of the battery
thd = 0; % threshold counter for U_MPP

delta_ip = -0.01;

SOC_ip = 1:delta_ip:0; % interpolated state of charge of the battery
U_0_D_ip = interp1(SOC, U_0_D, SOC_ip); % interpolated open-circuit voltage of the battery when discharging
U_0_C_ip = interp1(SOC, U_0_C, SOC_ip); % interpolated open-circuit voltage of the battery when charging
R_e_D_ip = interp1(SOC, R_e_D, SOC_ip); % interpolated electrolyte resistance of the battery when discharging
R_e_C_ip = interp1(SOC, R_e_C, SOC_ip); % interpolated electrolyte resistance of the battery when charging

U_0 = (U_0_D_ip + U_0_C_ip) / 2; % mean open circuit voltage

for n = 1:N_md % SCC limits power
    for o = 1:length(t_S)
        if(P_MPP(n,o) > P_MPP_max)
            P_MPP(n,o) = P_MPP_max;
            I_MPP(n,o) = P_MPP(n,o) / U_MPP(n,o);
        end
    end
end

for n = 1:N_md % PV voltage at SCC input terminal after cable losses
    for o = 1:length(t_S)
       U_in(n,o) = U_MPP(n,o) - 2 * (R_wires(1) + R_wires(2)) * I_MPP(n,o);
    end
end

[~,idx_SOC] = min(abs(SOC_ip - SOC_init)); % find closest index

for n = 1:N_md % calculating SOC and voltage of the battery (improve in version 1.1)
   for o = 1:length(t_S)
        if(n == 1 && o < idx_mission_se(n,1)) % mission on the first day has not started yet 
            thd = -1;
        elseif(n == 1 && o == idx_mission_se(n,1)) % if the mission day starts at 0:00
            if(U_in(n,o) > (U_0(idx_SOC) + 5) && (thd == 0 || thd == -1)) % I_MPP on
                thd = 1;
            elseif(U_in(n,o) < (U_0(idx_SOC) + 1) && thd == 1) % I_MPP off;
                thd = 0;
            end
        else
            if(o - 1 == 0) % if a new day started, get U_B from the last time entry from the previous day
                if(U_in(n,o) > (U_B(n - 1,length(t_S)) + 5) && (thd == 0 || thd == -1)) % I_MPP on
                    thd = 1;
                elseif(U_in(n,o) < (U_B(n - 1,length(t_S)) + 1) && thd == 1) % I_MPP off;
                    thd = 0;
                end
            else % get U_B from the same day but from the previous time entry
                if(U_in(n,o) > (U_B(n,o - 1) + 5) && (thd == 0 || thd == -1)) % I_MPP on
                    thd = 1;
                elseif(U_in(n,o) < (U_B(n,o - 1) + 1) && thd == 1) % I_MPP off;
                    thd = 0;
                end
            end 
        end
        
        if(thd == 1) % PV generator is connected (improve in version 1.1)
            I_B(n,o) =  I_SCC - I_MPP(n,o) * eta_SCC + I_L(n,o);
        elseif(thd == 0) % PV generator is disconnected
            I_B(n,o) =  I_SCC + I_L(n,o);
        elseif(thd == -1) % before the mission starts on the first day
            I_B(n,o) = 0;
        end
        
        if(I_B(n,o) > I_B_max)
            I_B(n,o) = I_B_max;
        elseif(I_B(n,o) < - I_B_max)
            I_B(n,o) = - I_B_max;
        end
        
        if(I_B(n,o) > 0) % battery is discharging 
            Q_tot = Q_nom * (I_nom / I_B(n,o))^(k_P - 1);
            
            if(n > 1 && o == 1) % get SOC from the end of the previous day
                SOC_curr(n,o) = SOC_curr(n - 1, length(t_S)) - (1 / Q_tot) * I_B(n,o) * delta_t_S;
            else % get previous state of charge
                SOC_curr(n,o) = SOC_curr(n,o - 1) - (1 / Q_tot) * I_B(n,o) * delta_t_S;
            end
            
            if(SOC_curr(n,o) > 1)
                SOC_curr(n,o) = 1;
            elseif(SOC_curr(n,o) < 0)
                SOC_curr(n,o) = 0;
            end
            
            [~,idx_SOC] = min(abs(SOC_ip - SOC_curr(n,o))); % find closest index
            c_B = abs(U_0_C_ip(idx_SOC) - U_0_D_ip(idx_SOC)) / 2;
            U_B(n,o) = U_0(idx_SOC) - c_B - R_e_D_ip(idx_SOC) * I_B(n,o);
            
        elseif(I_B(n,o) < 0) % battery is charging 
            Q_tot = Q_nom;
            
            if(n > 1 && o == 1) % get SOC from the end of the previous day
                SOC_curr(n,o) = SOC_curr(n - 1, length(t_S)) - 1 / Q_tot * I_B(n,o) * delta_t_S;
            else % get previous state of charge
                SOC_curr(n,o) = SOC_curr(n,o - 1) - 1 / Q_tot * I_B(n,o) * delta_t_S;
            end
            
            if(SOC_curr(n,o) > 1)
                SOC_curr(n,o) = 1;
            elseif(SOC_curr(n,o) < 0)
                SOC_curr(n,o) = 0;
            end
            
            [~,idx_SOC] = min(abs(SOC_ip - SOC_curr(n,o))); % find closest index
            c_B = abs(U_0_C_ip(idx_SOC) - U_0_D_ip(idx_SOC)) / 2;
            U_B(n,o) = U_0(idx_SOC) + c_B - R_e_C_ip(idx_SOC) * I_B(n,o);
            
        elseif(I_B(n,o) == 0) % battery is in standby (SOC_init stays the same)
            if(n == 1 && o == 1)
                SOC_curr(n,o) = SOC_init;
            elseif(n > 1 && o == 1) % get SOC from the end of the previous day
                SOC_curr(n,o) = SOC_curr(n - 1,length(t_S));
            else % get previous state of charge
                SOC_curr(n,o) = SOC_curr(n,o - 1);
            end
            
            if(SOC_curr(n,o) > 1)
                SOC_curr(n,o) = 1;
            elseif(SOC_curr(n,o) < 0)
                SOC_curr(n,o) = 0;
            end
            
            [~,idx_SOC] = min(abs(SOC_ip - SOC_curr(n,o))); % find closest index
            U_B(n,o) = U_0(idx_SOC);
        end 
        
        U_L(n,o) =  U_B(n,o) - 2 * R_wires(3) * I_B(n,o); % voltage at repeater after cable losses
        if(U_L(n,o) < U_L_min)
            error('Voltage at repeater radio system is too low!');
        elseif(U_L(n,o) > U_L_max)
            error('Voltage at repeater radio system is too high!');
        end
   end
end

if(min(min(SOC_curr)) == 0) % battery got fully discharged
    error('LiFePO_4 battery gets fully discharged! Use more PV generators in parallel!');
end

%% [OUTPUT] Command window and plots

fnct_cw_output(show_command_window_output, N_md, gamma_S, alpha_S, ...
    date_array, beta, t_S, W_G_el, SOC_curr); % command window output
 
diary off
