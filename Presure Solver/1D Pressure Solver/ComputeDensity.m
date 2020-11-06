
function [rho, drhodp] = ComputeDensity(p, po,rhoo,c_rho)

rho = rhoo.*(exp(c_rho .*(p -po)));
drhodp = c_rho.*rho;
end 
