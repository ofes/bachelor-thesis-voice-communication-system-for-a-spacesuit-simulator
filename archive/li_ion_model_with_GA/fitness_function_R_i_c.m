%% Charging R_i fitness function
function R_i = fitness_function_R_i_c(a)
SOC = 0.1;
C_r = 0.5;
R_i  = (a(1) + a(2) * C_r + a(3) * C_r^2) * exp(- a(4) * SOC) + (a(5) + a(6) * C_r + a(7) * C_r^2);
end