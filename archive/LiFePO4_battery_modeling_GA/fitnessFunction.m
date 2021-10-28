%% Fitness function for the battery terminal voltage
% t is in (s)
% Q_rem is in (As)

function y = fitnessFunction(a)

    I_D = 1.5;                              % Discharge current
    Q_rem = 3600 * [50 45 27.5 5 0];        % Remaining capacity in (As)
    DOD = [0 0.1 0.45 0.9 1];               % Depth of discharge
    U_B_M = [14.4 13.25 13.20 13.13 10];    % Measured voltages from data sheet
    t = 3600 * [0 3.33 15 30 33.3];         % Time in (s)
    U_B_C = zeros(size(U_B_M));
    R_1 = zeros(size(U_B_C));
    R_2 = zeros(size(U_B_C));
    C = zeros(size(U_B_C));
    U_0 = zeros(size(U_B_C));
    x = zeros(size(Q_rem));
    z = zeros(size(DOD));
    tmp = 0;

    for i = 1:length(Q_rem)
        x(i) = I_D/Q_rem(i);
    end

    for i = 1:length(DOD)
        z(i) = 1 - DOD(i);
    end

    for i = 1:length(U_B_M)
        R_1(i) = ( a(1) + a(2) * x(i) + a(3) * x(i)^2 ) * exp( -a(4) * z(i) ) + ( a(5) + a(6) * x(i) + a(7) * x(i)^2 );
        R_2(i) = ( a(8) + a(9) * x(i) + a(10) * x(i)^2 ) * exp( -a(11) * z(i) ) + ( a(12) + a(13) * x(i) + a(14) * x(i)^2 );
        C(i) = - ( a(15) + a(16) * x(i) + a(17) * x(i)^2 ) * exp( -a(18) * z(i) ) + ( a(19) + a(20) * x(i) + a(21) * x(i)^2 );
        U_0(i) = ( a(22) + a(23) * x(i) + a(24) * x(i)^2 ) * exp( -a(25) * z(i) ) + ( a(26) + a(27) * z(i) + a(28) * z(i)^2 + a(29) * z(i)^3) - a(30) * x(i) + a(31) * x(i)^2;
    end
    
    for i = length(U_B_M)
        U_B_C(i) = ( Q_rem(i) / C(i) + I_D * R_2(i) ) * exp( -t(i) / ( C(i) * R_2(i) ) ) + U_0(i) - I_D * ( R_1(i) + R_2(i) );
    end
    
    for i = 1:length(U_B_M)
        tmp = tmp + abs( U_B_M(i) - U_B_C(i) );
    end
    
    y = tmp^2;
end