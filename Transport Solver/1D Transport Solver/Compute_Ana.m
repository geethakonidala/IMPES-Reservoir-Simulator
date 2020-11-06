function [x_ana, sw2] = Compute_Ana(t,no,nw,mu_w,mu_o,sor,swc,krwe,kroe,phi)

ut = 1;
nsw = 100;                  % number of gridpoints
dsw = (1-sor-swc)/(nsw-1);
sw = swc:dsw:(1-sor);       % saturation range

% relperms + mobility + fracflow
krw = krwe*((sw-swc)/(1-swc-sor)).^nw;
kro = kroe*((1-sw-sor)/(1-swc-sor)).^no;
mob_w = krw/mu_w;
mob_o = kro/mu_o;
fw = mob_w./(mob_w+mob_o);

% fracflow derivative (numerical)
rr = 2:(nsw-1);
dfw1 = [0 (fw(rr+1)-fw(rr-1))./(sw(rr+1)-sw(rr-1)) 0];
vw = ut/phi*dfw1;

% jump velocity
dfw2 = (fw - fw(1))./(sw-sw(1));
[shock_vel, shock_index] = max(dfw2);
shock_sat = sw(shock_index);
rr2 = shock_index:nsw;

% valid saturation range within BL solution
sw2 = [sw(1) sw(1) sw(rr2)];
vw2 = vw(rr2);
vw2 = [1.5*vw2(1) vw2(1) vw2];

x_ana = vw2*t;


end 