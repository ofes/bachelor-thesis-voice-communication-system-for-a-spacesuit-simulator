%% Main: PV generator model
% This script plots the current-voltage and power-voltage characteristic of the PV generator used to supply the self-sufficient repeater radio infrastructure, depending on the solar irradiance and the PV cell temperature. 

% Organization:     OeWF (Austrian Space Forum)
% Author:           Omar Filip El Sendiouny
% Project:          Serenity BU-COMMs
% Date:             10.04.2021
% Version:          1.0

clear all;
close all;
clc;

figure_counter = 1;
fnt_sz = 17.5;

%% PV generator specifications:
% MPP ... maximim power point
% STC ... standard test conditions
% SC ... short-circuit
% OC ... open-circuit
% Input style: [PV_generator_1 PV_generator_2  ...]

number_PV_generators = 2; % [AE_Solar_AE195SMM6-36 DAS_Energy_DAS145PF]

N_C = [36 36]; % number of PV cells in (1)
I_MPP_STC = [9.11 8.19]; % current at MPP for STC in (A)
U_MPP_STC = [21.41 18.69]; % voltage at MPP for STC in (V)
I_SC_STC = [9.79 8.69]; % SC current for STC in (A)
U_OC_STC = [24.27 22.85]; % OC voltage for STC in (V)
TC_I_SC = [0.05 0.051]; % temperature coefficient for I_SC in (%/degrees C)
TC_U_OC = [-0.29 -0.31]; % temperature coefficient for U_OC in (%/degrees C)
m = [1.19045 1.58972]; % ideality factor
pv_name = {'ae_solar', 'das_energy'}; % for file exporting

E_STC = 1000; % solar irradiance received by an inclined PV generator for STC in (W/m(h)^2)
theta_STC = 25; % ambient temperature for STC in (degrees C)
e = 1.602176634 * 10^(-19); % elementary charge in (As)
k_B = 1.380649 * 10^(-23); % Bolzmann constant in (Ws/K)

for h = 1:number_PV_generators
    %% Plotting the effect of the solar irradiance onto an inclined PV generator:
    % In total seven different irradiance values were examined.

    theta_C = [25 25 25 25 25 25]; % PV cell temperature for the PV generator in (degrees C)
    E_G = [200 400 600 800 1000 1200]; % solar irradiance received by an inclined PV generator in (W/m(h)^2)

    % ----------------------------------------------------------------- %
    % Current-voltage characteristic                                    %
    % ----------------------------------------------------------------- %

    res = 0.0001; % resolution

    for i = length(E_G):-1:1

        U_T = k_B * (theta_C(i) + 273.15)/e; % calc. U_T(theta_C) (V)
        I_Ph = I_SC_STC(h) * E_G(i) / E_STC * (1 + TC_I_SC(h) / 100 * (theta_C(i) - theta_STC)); % calc. I_Ph(theta_C, E_G) in (A)
        I_Ph_const_irr = I_SC_STC(h) * (1 + TC_I_SC(h) / 100 * (theta_C(i) - theta_STC)); % calc. I_Ph(theta_C, E_STC) in (A)

        U_OC_0 = U_OC_STC(h) * (1 + TC_U_OC(h) / 100 * (theta_C(i) - theta_STC)) + m(h) * N_C(h) * U_T * log(I_Ph/I_Ph_const_irr); % calc. U_OC(theta_C, E_G) in (V)

        I_S_0 = I_Ph * exp(-U_OC_0/(m(h) * N_C(h) * U_T)); % calc. I_S(theta_C) in (A)

        U_PV = 0:res:U_OC_0; % init. U_PV in (V)
        I_PV = zeros(1, length(U_PV)); % init. I_PV in (A)

        for j = 1:length(U_PV)
            I_PV(j) =  I_Ph + I_S_0 - I_S_0 * exp(U_PV(j)/(m(h) * N_C(h) * U_T)); % calc. I_PV(U_PV, theta_C, Phi_G)
        end 

        figure(figure_counter); % plotting results
        grid on;
        hold all;
        disp_lgnd = ['$E_\mathrm{G} =$ ', num2str(E_G(i)),'$\mathrm{Wm^{-2}}$'];
        plot(U_PV, I_PV, 'DisplayName', disp_lgnd, 'Linewidth', 1.4);

        xlabel('PV generator voltage $U_\mathrm{PV}$ in $\left(\mathrm{V}\right)$', 'Interpreter', 'latex', 'FontSize', fnt_sz);
        ylabel('PV generator current $I_\mathrm{PV}$ in $\left(\mathrm{A}\right)$', 'Interpreter', 'latex', 'FontSize', fnt_sz);

        legend('-DynamicLegend', 'Location', 'SouthOutside', 'Interpreter', 'latex', 'FontSize', fnt_sz, 'NumColumns',3);
        ax = gca;
        set(ax,'TickLabelInterpreter', 'latex', 'FontSize', fnt_sz);

        plot_pos_x0 = 10;
        plot_pos_y0 = 10;
        plot_width = 600;
        plot_height = 450;
        set(gcf,'position',[plot_pos_x0, plot_pos_y0, plot_width, plot_height]);
    end

    export_title = ['latex_export/image_curr_volt_irr_',  pv_name{h},'.eps'];
    exportgraphics(gcf, export_title);
    figure_counter = figure_counter + 1;

    % ----------------------------------------------------------------- %
    % Power-voltage characteristic                                      %
    % ----------------------------------------------------------------- %

    for i = length(E_G):-1:1

        U_T = k_B * (theta_C(i) + 273.15)/e;
        I_Ph = I_SC_STC(h) * E_G(i) / E_STC * (1 + TC_I_SC(h) / 100 * (theta_C(i) - theta_STC));
        I_Ph_const_irr = I_SC_STC(h) * (1 + TC_I_SC(h) / 100 * (theta_C(i) - theta_STC));

        U_OC_0 = U_OC_STC(h) * (1 + TC_U_OC(h) / 100 * (theta_C(i) - theta_STC)) + m(h) * N_C(h) * U_T * log(I_Ph/I_Ph_const_irr);

        I_S_0 = I_Ph * exp(-U_OC_0/(m(h) * N_C(h) * U_T));

        U_PV = 0:res:U_OC_0;
        P_PV = zeros(1, length(U_PV)); % init. P_PV in (W)

        for j = 1:length(U_PV)
            P_PV(j) = U_PV(j) * (I_Ph + I_S_0 - I_S_0 * exp(U_PV(j)/(m(h) * N_C(h) * U_T))); % calc. P_PV(U_PV, theta_C, Phi_G)
        end 

        figure(figure_counter); % plotting results
        grid on;
        hold all;
        disp_lgnd = ['$E_\mathrm{G} =$ ', num2str(E_G(i)),'$\mathrm{Wm^{-2}}$'];
        plot(U_PV, P_PV, 'DisplayName', disp_lgnd, 'Linewidth', 1.4);

        xlabel('PV generator voltage $U_\mathrm{PV}$ in $\left(\mathrm{V}\right)$', 'Interpreter', 'latex', 'FontSize', fnt_sz);
        ylabel('PV generator power $P_\mathrm{PV}$ in $\left(\mathrm{W}\right)$', 'Interpreter', 'latex', 'FontSize', fnt_sz);

        legend('-DynamicLegend', 'Location', 'SouthOutside', 'Interpreter', 'latex', 'FontSize', fnt_sz, 'NumColumns',3);
        ax = gca;
        set(ax,'TickLabelInterpreter', 'latex', 'FontSize', fnt_sz);

        plot_pos_x0 = 10;
        plot_pos_y0 = 10;
        plot_width = 600;
        plot_height = 450;
        set(gcf,'position',[plot_pos_x0, plot_pos_y0, plot_width, plot_height]);
    end

    export_title = ['latex_export/image_power_volt_irr_',  pv_name{h},'.eps'];
    exportgraphics(gcf, export_title);
    figure_counter = figure_counter + 1;

    %% Plotting the effect of the solar irradiance onto an inclined PV generator:
    % In total seven different irradiance values were examined.

    theta_C = [-50 -25 0 25 50 75]; % PV cell temperature for the PV generator in (degrees C)
    E_G = [1000 1000 1000 1000 1000 1000]; % solar irradiance received by an inclined PV generator in (W/m(h)^2)

    % ----------------------------------------------------------------- %
    % Current-voltage characteristic                                    %
    % ----------------------------------------------------------------- %

    res = 0.0001;

    for i = length(E_G):-1:1

        U_T = k_B * (theta_C(i) + 273.15)/e;
        I_Ph = I_SC_STC(h) * E_G(i) / E_STC * (1 + TC_I_SC(h) / 100 * (theta_C(i) - theta_STC));
        I_Ph_const_irr = I_SC_STC(h) * (1 + TC_I_SC(h) / 100 * (theta_C(i) - theta_STC));

        U_OC_0 = U_OC_STC(h) * (1 + TC_U_OC(h) / 100 * (theta_C(i) - theta_STC)) + m(h) * N_C(h) * U_T * log(I_Ph/I_Ph_const_irr); 

        I_S_0 = I_Ph * exp(-U_OC_0/(m(h) * N_C(h) * U_T));

        U_PV = 0:res:U_OC_0;
        I_PV = zeros(1, length(U_PV));

        for j = 1:length(U_PV)
            I_PV(j) =  I_Ph + I_S_0 - I_S_0 * exp(U_PV(j)/(m(h) * N_C(h) * U_T));
        end 

        figure(figure_counter); % plotting results
        grid on;
        hold all;
        disp_lgnd = ['$\vartheta_\mathrm{C} =$ ', num2str(theta_C(i)),'$^\circ\mathrm{C}$'];
        plot(U_PV, I_PV, 'DisplayName', disp_lgnd, 'Linewidth', 1.4);

        xlabel('PV generator voltage $U_\mathrm{PV}$ in $\left(\mathrm{V}\right)$', 'Interpreter', 'latex', 'FontSize', fnt_sz);
        ylabel('PV generator current $I_\mathrm{PV}$ in $\left(\mathrm{A}\right)$', 'Interpreter', 'latex', 'FontSize', fnt_sz);

        legend('-DynamicLegend', 'Location', 'SouthOutside', 'Interpreter', 'latex', 'FontSize', fnt_sz, 'NumColumns',3);
        ax = gca;
        set(ax,'TickLabelInterpreter', 'latex', 'FontSize', fnt_sz);

        plot_pos_x0 = 10;
        plot_pos_y0 = 10;
        plot_width = 600;
        plot_height = 450;
        set(gcf,'position',[plot_pos_x0, plot_pos_y0, plot_width, plot_height]);
    end

    export_title = ['latex_export/image_curr_volt_temp_',  pv_name{h},'.eps'];
    exportgraphics(gcf, export_title);
    figure_counter = figure_counter + 1;

    % ----------------------------------------------------------------- %
    % Power-voltage characteristic                                      %
    % ----------------------------------------------------------------- %

    for i = length(E_G):-1:1

        U_T = k_B * (theta_C(i) + 273.15)/e;
        I_Ph = I_SC_STC(h) * E_G(i) / E_STC * (1 + TC_I_SC(h) / 100 * (theta_C(i) - theta_STC));
        I_Ph_const_irr = I_SC_STC(h) * (1 + TC_I_SC(h) / 100 * (theta_C(i) - theta_STC));

        U_OC_0 = U_OC_STC(h) * (1 + TC_U_OC(h) / 100 * (theta_C(i) - theta_STC)) + m(h) * N_C(h) * U_T * log(I_Ph/I_Ph_const_irr);

        I_S_0 = I_Ph * exp(-U_OC_0/(m(h) * N_C(h) * U_T));

        U_PV = 0:res:U_OC_0;
        P_PV = zeros(1, length(U_PV));

        for j = 1:length(U_PV)
            P_PV(j) = U_PV(j) * (I_Ph + I_S_0 - I_S_0 * exp(U_PV(j)/(m(h) * N_C(h) * U_T)));
        end 
        
        if theta_C(i) == 25 % used to adjust the ideality factor
            disp_str = ['PV generator: ', pv_name{h}];
            disp(disp_str)
            disp_str = ['PV cell temperature = ', num2str(theta_C(i)), 'degrees C.'];
            disp(disp_str)
            disp('E_G = 1000Wm^(-2).')
            disp_str = ['Output power at MPP = ', num2str(max(P_PV)), 'W.'];
            disp(disp_str)
            disp_str = ['Ideality factor = ', num2str(m(h))];
            disp(disp_str)
            disp(' ')
        end

        figure(figure_counter); % plotting results
        grid on;
        hold all;
        disp_lgnd = ['$\vartheta_\mathrm{C} =$ ', num2str(theta_C(i)),'$^\circ\mathrm{C}$'];
        plot(U_PV, P_PV, 'DisplayName', disp_lgnd, 'Linewidth', 1.4);

        xlabel('PV generator voltage $U_\mathrm{PV}$ in $\left(\mathrm{V}\right)$', 'Interpreter', 'latex', 'FontSize', fnt_sz);
        ylabel('PV generator power $P_\mathrm{PV}$ in $\left(\mathrm{W}\right)$', 'Interpreter', 'latex', 'FontSize', fnt_sz);

        legend('-DynamicLegend', 'Location', 'SouthOutside', 'Interpreter', 'latex', 'FontSize', fnt_sz, 'NumColumns',3);
        ax = gca;
        set(ax,'TickLabelInterpreter', 'latex', 'FontSize', fnt_sz);

        plot_pos_x0 = 10;
        plot_pos_y0 = 10;
        plot_width = 600;
        plot_height = 450;
        set(gcf,'position',[plot_pos_x0, plot_pos_y0, plot_width, plot_height]);
    end

    export_title = ['latex_export/image_power_volt_temp_',  pv_name{h},'.eps'];
    exportgraphics(gcf, export_title);
    figure_counter = figure_counter + 1;
end