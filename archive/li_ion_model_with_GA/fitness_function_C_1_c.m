%% Charging R_1 fitness function
function C_1 = fitness_function_C_1_c(a)
SOC = 0.5;
C_r = 1;
C_1  = - (a(15) + a(16) * C_r + a(17) * C_r^2) * exp(- a(18) * SOC) + (a(19) + a(20) * C_r + a(21) * C_r^2);
end