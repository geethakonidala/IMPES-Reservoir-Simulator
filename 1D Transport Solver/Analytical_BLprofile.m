%----------------------------------------------------
% ANALITYCAL SOLUTION OF BUCKLEY LEVERETT EQUATION
% ----------------------------------------------------
clear all; close all;
% parameters
mu_w = 1.0e-3;
mu_o = 10.0e-3;
nw = 3.5;
no = 2;
ut = 1;
swc = 0.2;
sor = 0.1;
krwe = 0.7;
kroe = 0.8;
porosity= 0.3;

% saturation
nsw = 101; % number of gridpoints
dsw = (1-sor-swc)/(nsw-1);
sw = swc:dsw:(1-sor); % saturation range

% relperms + mobility + fracflow
krw = krwe*((sw-swc)/(1-swc-sor)).^nw;
kro = kroe*((1-sw-sor)/(1-swc-sor)).^no;
mob_w = krw/mu_w;
mob_o = kro/mu_o;
fw = mob_w./(mob_w+mob_o);

% fracflow derivative (numerical)
rr = 2:(nsw-1);
dfw1 = [0 (fw(rr+1)-fw(rr-1))./(sw(rr+1)-sw(rr-1)) 0];
vw = ut/porosity*dfw1;

% jump velocity
dfw2 = (fw - fw(1))./(sw-sw(1));
[shock_vel, shock_index] = max(dfw2);
shock_sat = sw(shock_index);
rr2 = shock_index:nsw;

% valid saturation range within BL solution
sw2 = [sw(1) sw(1) sw(rr2)];
vw2 = vw(rr2);
vw2 = [1.5*vw2(1) vw2(1) vw2];
t = 10;
x = vw2*t;

% plot of saturation as function of position
figure(1);
subplot(2,1,2); plot(x,sw2,'.-','linewidth',2);
xlabel('x [m]'); ylabel('S_w [-]');
% plot of fractional flow derivative + jump velocity
subplot(2,1,1); plot(sw,dfw1,sw,dfw2,[shock_sat shock_sat],[0
max(dfw1)],'r:',[swc (1-sor)],[shock_vel shock_vel],'r:');
xlabel('Sw'); ylabel('derivative');
legend('d{f_w}/d{S_w}','\Deltaf_w/\DeltaS_w');