function [u_vosst, v_vosst] = vosstanovitelnaya_function( k, dk, h, l, veter0)


opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-12, 'MaxStep', (h(end)-h(1))/1000);

[~, y] = ode15s(@(t, y) odefun(t, y), h, veter0, opt);
veter0

u_vosst = y(:,1)';
v_vosst = y(:,3)';

function dydx = odefun(x,y)
    
kz = ppval(k, x);
dkz = ppval(dk, x);


dydx = [y(2);
    -(y(2)*dkz + l*y(3)  ) / kz;
    y(4);
    -(y(4)*dkz - l*y(1)  ) / kz];

end

end