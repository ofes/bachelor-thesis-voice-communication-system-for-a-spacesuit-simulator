figure(1)
plot(t_S, theta)
figure(2)
plot(t_S, gamma_S)
figure(3)
plot(t_S, alpha_S)


r = cos(gamma_S * pi /180);
figure(4)
polarplot(deg2rad(alpha_S),r,'Linewidth', 1.4, 'Color', '#f5a742');
rlim([0 1]);
pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'left';
thetaticks([0 90 180 270])
thetaticklabels({'South: \alpha_S = 0°','West: \alpha_S = 90°','North: \alpha_S= 180°','East: \alpha_S = -90°'})
rticklabels({})
rticks([])