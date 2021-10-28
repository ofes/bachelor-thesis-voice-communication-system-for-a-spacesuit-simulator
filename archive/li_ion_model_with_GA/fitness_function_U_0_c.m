%% Charging R_1 fitness function
function U_0 = fitness_function_U_0_c(a)SOC = 0.5;
C_r = 1;
U_0  = (a(22) + a(23) * C_r + a(24) * C_r^2) * exp(- a(25) * SOC) + (a(26) + a(27) * SOC + a(28) * SOC^2 + a(29) * SOC^3) - a(30) * C_r + a(31) * C_r^2;
end