%% Main: Attenuation curve
% This script plots the attenuation at 20 degrees C for the Messi & Paoloni Ultrafelx 7 coaxial cable. 

% Organization:     OeWF (Austrian Space Forum)
% Author:           Omar Filip El Sendiouny
% Project:          Serenity BU-COMMs
% Date:             05.04.2021
% Version:          1.0

clear all;
close all;
clc;

fnt_sz = 17.5;

%% Datasheet extraction and interpolation

f = [1.8 3.5 7 10 14 21 28 50 100 144 200 400 430 800 1000 1296 2400 3000 4000 5000 6000 7000 8000]; % frequency in (MHz) extracted from datasheet
L_UF7 = [1.1 1.3 1.7 1.9 2.2 2.6 3.0 4.0 5.8 6.9 8.2 11.8 12.3 17.1 19.3 22.3 32.3 36.2 42.6 49.3 55.3 61.6 68.4]; % attenuation in (dB/100m) extracted from datasheet

f_int = 2:0.05:8000; % generating frequency vector (0.05MHz resolution)
L_UF7_int = interp1(f, L_UF7, f_int); % interpolation of attenuation data

%% Obtaining attenuation for 158,950MHz (current OeWF VHF frequency)

f_target = 158.950; % target frequency for which the attenuation needs to be found

for i = 1:length(f_int) % find target frequency
    if(round(f_int(i), 3) == f_target) % rounding required necessary due to initialization problems with f_int
        print_str = ['The attenuation for f = ',num2str(f_target),'MHz is ', num2str(L_UF7_int(i)), 'dB/100m.'];
        disp(print_str);
        break;
    end
end

%% Plotting the interpolated curve

figure(1);
grid on;
hold on;
plot(f_int, L_UF7_int, 'Linewidth', 1.4, 'Color', '#a11b1b')

xlabel('Frequency $f$ in $\left(\mathrm{MHz}\right)$', 'Interpreter', 'latex', 'FontSize', fnt_sz);
ylabel('Attenuation $L_\mathrm{UF7}$ in $\left(\mathrm{dB}/100\mathrm{m}\right)$', 'Interpreter', 'latex', 'FontSize', fnt_sz);

ax = gca;
set(ax,'TickLabelInterpreter', 'latex', 'FontSize', fnt_sz, 'XScale', 'log');

plot_pos_x0 = 10;
plot_pos_y0 = 10;
plot_width = 600;
plot_height = 400;
set(gcf,'position',[plot_pos_x0, plot_pos_y0, plot_width, plot_height]);

export_title = 'latex_export/image_attenuation_UF7.eps';
exportgraphics(gcf, export_title);

