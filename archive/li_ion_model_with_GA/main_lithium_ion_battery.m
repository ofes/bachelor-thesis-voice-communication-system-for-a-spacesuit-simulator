%% LiFePO4 battery modeling:
% The EEC of the battery is modeled with the first-order Thevenin model.
% The missiong lumped-parameters are are obtained by useing the Genetic
% Algorithm.

clc;
clear all;
close all;

Funct_R_i = @fitness_function_R_i_c;
Funct_R_1 = @fitness_function_R_1_c;
Funct_C_1 = @fitness_function_C_1_c;
Funct_U_0 = @fitness_function_U_0_c;

nvars_R_i = 7;
nvars_R_1 = 7;
nvars_C_i = 7;
nvars_U_0 = 10;

[a, R_i] = ga(Funct_R_i, nvars_R_i);

a
R_i

[b, R_1] = ga(Funct_R_1, nvars_R_1);

b
R_1
