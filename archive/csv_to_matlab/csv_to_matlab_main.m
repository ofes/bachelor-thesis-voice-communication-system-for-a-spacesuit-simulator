clear all;
close all;
clc;

A = table2array(readtable('i_exp_dec.csv'));

sz = size(A);

x = zeros(1, sz(1));
y = zeros(1, sz(1));

% B = smoothdata(A);
% 
% for i = 1:2000
%     x(i) = B(i,1);
%     y(i) = B(i,2) * (-1);
% end

for i = 1:2000
    u(i) = A(i,1);
    v(i) = smoothdata(A(i,2) * (-1));
end

figure(1);
a = plot(u, v, 'Linewidth', 1.5);
hold on

% MC = ischange(w, 'variance', 'Threshold', 5000);
% count = 0;
% for i = 1:2000
%     if(MC(i)) 
%         count = count + 1;
%     end
% end
% 
% idxmax = find(MC == max(MC));
% 
% figure(2);
% yyaxis right;
% c = plot(u, w, 'Linewidth', 2);
% yyaxis left;
% %d = plot(u, MC, '*', 'MarkerIndices',[idxmax]);
% e = stem(u,MC);
% set(e, 'Marker', 'none')