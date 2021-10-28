%% Determining discharging curve with GA

close all;
clear all;
clc;

%% Finding coefficients

FitFcn = @fitnessFunction;
nvars = 31;
options.CrossoverFraction = 0.85;

[a, fval] = ga(FitFcn, nvars, options);

a
fval

%% Plotting discharge curve

I_D = 1.5;
delta_t = 1;
t = 0:delta_t:(3600 * 33.3);
Q_rem = zeros(size(t));
DOD = zeros(size(Q_rem));
x = zeros(size(Q_rem));
z = zeros(size(Q_rem));
U_B = zeros(size(Q_rem));
R_1 = zeros(size(U_B));
R_2 = zeros(size(U_B));
C = zeros(size(U_B));
U_0 = zeros(size(U_B));
Q_0 = 3600 * 50;

Q_rem(1) = Q_0;
for i = 2:1:length(t)
    Q_rem(i) = Q_rem(i-1) - I_D * delta_t;
end

for i = 1:length(Q_rem)
    DOD(i) = Q_rem(i)/Q_0;
end

for i = 1:length(Q_rem)
    z(i) = 1 - DOD(i);
end

for i = 1:length(Q_rem)
    x(i) = I_D/Q_rem(i);
end

for i = 1:length(Q_rem)
    R_1(i) = ( a(1) + a(2) * x(i) + a(3) * x(i)^2 ) * exp( -a(4) * z(i) ) + ( a(5) + a(6) * x(i) + a(7) * x(i)^2 );
    R_2(i) = ( a(8) + a(9) * x(i) + a(10) * x(i)^2 ) * exp( -a(11) * z(i) ) + ( a(12) + a(13) * x(i) + a(14) * x(i)^2 );
    C(i) = - ( a(15) + a(16) * x(i) + a(17) * x(i)^2 ) * exp( -a(18) * z(i) ) + ( a(19) + a(20) * x(i) + a(21) * x(i)^2 );
    U_0(i) = ( a(22) + a(23) * x(i) + a(24) * x(i)^2 ) * exp( -a(25) * z(i) ) + ( a(26) + a(27) * z(i) + a(28) * z(i)^2 + a(29) * z(i)^3) - a(30) * x(i) + a(31) * x(i)^2;
end

for i = 1:length(Q_rem)
    U_B(i) = ( Q_rem(i) / C(i) + I_D * R_2(i) ) * exp( -t(i) / ( C(i) * R_2(i) ) ) + U_0(i) - I_D * ( R_1(i) + R_2(i) );
end

plot(DOD, U_B);












