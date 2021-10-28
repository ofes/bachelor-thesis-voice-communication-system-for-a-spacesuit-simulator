clear all;
close all;
clc;

CF = 2 * pi / 360;

delta = -23.45 * CF;
phi = 48.2 * CF;
hs(2) = 0;

h(3) = acos(-tan(delta) * tan(phi));
h(1) = - acos(-tan(delta) * tan(phi));

gamma_s = asin(sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(h));

for i = 1:3
    if h(i) < 0 
        alpha_s(i) = - acos( ( sin(phi) * cos(delta) * cos(h(i)) - cos(phi) * sin(delta) ) / ( cos(gamma_s(i)) ) );
    else
        alpha_s(i) = acos( ( sin(phi) * cos(delta) * cos(h(i)) - cos(phi) * sin(delta) ) / ( cos(gamma_s(i)) ) );
    end
end

fprintf('Solar sunrise hour angle: \t h_s,r = %g°\n', round(h(1)*1/CF, 2));
fprintf('Solar noon hour angle: \t\t h_s,n = %g°\n', round(h(2)*1/CF, 2));
fprintf('Solar sunset hour angle: \t h_s,r = %g°\n', round(h(3)*1/CF, 2));
fprintf('\n');

fprintf('Solar sunrise altitude: \t gamma_s,r = %g°\n', round(gamma_s(1)*1/CF, 2));
fprintf('Solar noon altitude: \t\t gamma_s,n = %g°\n', round(gamma_s(2)*1/CF, 2));
fprintf('Solar sunset altitude: \t\t gamma_s,s = %g°\n', round(gamma_s(3)*1/CF, 2));
fprintf('\n');

fprintf('Solar sunrise azimuth: \t\t alpha_s,r = %g°\n', round(alpha_s(1)*1/CF, 2));
fprintf('Solar noon azimuth: \t\t alpha_s,n = %g°\n', round(alpha_s(2)*1/CF, 2));
fprintf('Solar sunset azimuth: \t\t alpha_s,s = %g°\n', round(alpha_s(3)*1/CF, 2));


