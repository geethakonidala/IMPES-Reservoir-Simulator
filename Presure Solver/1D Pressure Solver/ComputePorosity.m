
function [phi, dphidp] = ComputePorosity(p, po,phio,ceff)
phi = phio.*(exp(ceff .*(p -po)));
dphidp = ceff.*phi;
end 