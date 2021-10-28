%% function_PVG_alignment

% Computes how a photovoltaic generator needs to be aligned at a given
% location on Earth. 

% Organization:     OeWF (Austrian Space Forum)
% Author:           Omar Filip El Sendiouny
% Project:          Serenity BU-COMMs
% Date:             27.12.2020
% Version:          0.1

% INPUT: lat = latitude
% INPUT: lon = longitude
% INPUT: date (vector: mission start, mission end)

% OUTPUT: beta = PVG's angle of inclination for each mission day (vector).
% OUTPUT: alpha_S = Sun's azimuth angle for each mission day (vector).
% OUTPUT: lat_Z = Sun's zenith latitude for each mission day (vector).
% OUTPUT: delta = Sun's declination for each mission day (vector).

function [beta, alpha_S, lat_Z, delta] = function_PVG_alignment (lat, lon, date)
    
    n_d = zeros(1, days(date(2) - date(1)) + 1);        % Initializing vector. Number of days since January 1st.
    delta = zeros(1, days(date(2) - date(1)) + 1);      % Initializing vector.
    gamma_S = zeros(1, days(date(2) - date(1)) + 1);    % Initializing vector. Sun's altitude.
    alpha_S = zeros(1, days(date(2) - date(1)) + 1);    % Initializing vector.
    z_h = zeros(1, days(date(2) - date(1)) + 1);        % Initializing vector. Equation of time.
    
    cf = 2 * pi / 360;                                  % Correction factor (deg to rad).
    current_date = date(1);                             % Mission start.
    
    for i = 1:(days(date(2) - date(1)) + 1)
        n_d(i) = 30.3 * (month(current_date) - 1) + day(current_date);
        delta(i) = 23.45 * sin(360 * (284 + n_d(i))/(365) * cf);
        current_date = datetime(current_date + days(1), 'Format', 'd-MMM-y');
    end
   
    % FINISH
    
end