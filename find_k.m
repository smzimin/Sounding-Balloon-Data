function k = find_k( u, v, h, l )
du1 = fnder(u, 1);
du2 = fnder(u, 2);
du3 = fnder(u, 3);

dv1 = fnder(v, 1);
dv2 = fnder(v, 2);
dv3 = fnder(v, 3);

solinit = bvpinit(h, [1 0]);
sol = bvp4c(@odefun, @twobc, solinit);
k = deval(sol, h);


function res = twobc(ya,yb)

ua = ppval(u, h(1));
d1ua = ppval(du1, h(1));
d2ua = ppval(du2, h(1));

va = ppval(v, h(1));
d1va = ppval(dv1, h(1));
d2va = ppval(dv2, h(1));


ub = ppval(u, h(end));
d1ub = ppval(du1, h(end));
d2ub = ppval(du2, h(end));

vb = ppval(v, h(end));
d1vb = ppval(dv1, h(end));
d2vb = ppval(dv2, h(end));

res = [(d1ua^2+d1va^2)*ya(2) + (d2ua*d1ua+d2va*d1va)*ya(1) + l*(ua*d1va-va*d1ua) ;
    (d1ub^2+d1vb^2)*yb(2) + (d2ub*d1ub+d2vb*d1vb)*yb(1) + l*(ub*d1vb-vb*d1ub)];

end

function dydx = odefun(x,y)
    
d1uz = ppval(du1, x);
d2uz = ppval(du2, x);
d3uz = ppval(du3, x);

d1vz = ppval(dv1, x);
d2vz = ppval(dv2, x);
d3vz = ppval(dv3, x);

dydx = [y(2);
    -(2*(d2uz*d1uz+d2vz*d1vz)*y(2) + (d3uz*d1uz+d3vz*d1vz)*y(1) ) / (d1uz^2+d1vz^2)];

end

end