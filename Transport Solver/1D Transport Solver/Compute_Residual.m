function Residual = Compute_Residual(q, phi, dx,dt, f_w, Sw,Sw_new, u, N)
Residual= zeros(N,1);
%Residual of 1st cell
  Residual(1) = q(1) - (phi./dt).*(Sw_new(1) - Sw(1)) - (1/dx).* (f_w(1).* u(2));
 % since velocity is zero at interface,  so no fractional flow term to
 % the left side of the 1st cell
 
% Residual of remaining cells 
for i =2:N
Residual(i) = q(i).*f_w(i) - (phi/dt).*(Sw_new(i) - Sw(i)) - (1/dx).* (f_w(i).* u(i+1) - f_w(i-1).*u(i));
end 

end 


