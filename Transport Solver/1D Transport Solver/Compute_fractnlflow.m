function [f_w, dfw_ds ] = Compute_fractnlflow(mu_w, mu_o, Sw,krwe,kroe,sor,swc, nw, no, N)

% Corey Relative permeabilities
Krw   = krwe*((Sw-swc)/(1-swc-sor)).^nw;
Kro   = kroe*((1-Sw-sor)/(1-swc-sor)).^no;

%% Calculating Mobility & Fractional flow
Lamda_W     = Krw./mu_w; 
Lamda_O     = Kro./mu_o;
f_w         = Lamda_W./(Lamda_O + Lamda_W);

%% Calculating fractional flow - saturation derivative 
dLamda_W_ds = (nw*krwe./mu_w).*((Sw-swc)/(1-swc-sor)).^(nw-1)*(1/(1-swc-sor));
dLamda_O_ds = (no*kroe./mu_o).*((1-Sw-sor)/(1-swc-sor)).^(no-1)*(-1/(1-swc-sor));

dfw_ds = (dLamda_W_ds.*( Lamda_W + Lamda_O) - (dLamda_W_ds + dLamda_O_ds).*Lamda_W)./(Lamda_W + Lamda_O).^2; 

end


