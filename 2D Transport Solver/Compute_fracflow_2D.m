function [f_w, dfw_ds] = Compute_fracflow_2D( Sw, swc,krwe,kroe, sor, mu_w, mu_o,nw,no, Nx, Ny)


for i = 1:Ny
    for j = 1:Nx
            
Krw(i,j)        = krwe*((Sw(i,j) -swc)/(1-swc-sor)).^nw;
Kro(i,j)        = kroe*((1-Sw(i,j)-sor)/(1-swc-sor)).^no;

Lamda_W(i,j)   = Krw(i,j)./mu_w;
Lamda_O(i,j)     = Kro(i,j)./mu_o;

f_w = Lamda_W./(Lamda_O + Lamda_W);

%% Calculating fractional flow - saturation derivative 
dLamda_W_ds(i,j) = (nw*krwe./mu_w).*((Sw(i,j)-swc)/(1-swc-sor)).^(nw-1)*(1/(1-swc-sor));
dLamda_O_ds(i,j) = (no*kroe./mu_o).*((1-Sw(i,j)-sor)/(1-swc-sor)).^(no-1)*(-1/(1-swc-sor));

dfw_ds(i,j) = (dLamda_W_ds(i,j).*( Lamda_W(i,j) + Lamda_O(i,j)) - (dLamda_W_ds(i,j) + dLamda_O_ds(i,j)).*Lamda_W(i,j))./(Lamda_W(i,j) + Lamda_O(i,j)).^2; 
   
    end 
end
end 