%% LiFePO4 battery modeling

clc;
close all;
clear all;

Q_0 = 50;
U_0 = 12.8;
I_0 = 20;
k = 1.05;

C_1 = 50;
R_1 = 0.001;
R_e = 1;

I = 20;    % Entladestrom

Q = Q_0 * (I_0/I)^(k - 1);

t = 0:0.01:1000;
U = zeros(size(t));

for i = 1:1:length(t)
    U(i) = U_0 + Q/C_1 * exp( -t(i) / (R_1 * C_1) ) - I * R_e - I * R_1 * (1 - exp( -t(i) / (R_1 * C_1) ));
end

plot(t,U);