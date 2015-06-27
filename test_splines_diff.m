x = 1:10;
h = 1:0.5:10;
y = rand(1,10);

pp = spline(x,y);
p_der1 = fnder(pp, 1);
p_der2 = fnder(pp, 2);
p_der3 = fnder(pp, 3);
y_prime = ppval(pp, h);
dy_prime = ppval(p_der1, h);
d2y_prime = ppval(p_der2, h);
d3y_prime = ppval(p_der3, h);

plot(x, y, 'ob', h, y_prime, 'g', h, dy_prime, 'm', h, d2y_prime, 'r', h, d3y_prime, 'c');
% 
% s = interp1(x,y,h,'linear','extrap');
% plot(x, y, 'o', h, s, 'g'); hold on;
% s = interp1(x, y, h, 'spline','extrap' );
% plot(x, y, 'ob', h, s, 'm')

