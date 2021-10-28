%% Charging R_1 fitness function
function R_1 = fitness_function_R_1_c(b)
SOC = 0.1;
C_r = 0.5;
R_1  = (b(1) + b(2) * C_r + b(3) * C_r^2) * exp(- b(4) * SOC) + (b(5) + b(6) * C_r + b(7) * C_r^2);
end