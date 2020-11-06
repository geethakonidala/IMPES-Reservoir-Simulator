function R = Compute_Residual2D(q, phi, dx,dy,dt, f_w, Sw_n,Sw_nu, Vx,Vy, Nx,Ny)

%Residual of 1st in x axis and its corresponding row
  R(1,1) = q(1,1) - (phi./dt).*(Sw_nu(1,1) - Sw_n(1,1)) - (1/dx).* (f_w(1,1).* Vx(1,2)) -(1/dy).*(f_w(1,1).* Vy(2,1)); % E,N
 R(2:Ny,1) =  - (phi./dt).*(Sw_nu(2:Ny,1) - Sw_n(2:Ny,1)) - (1/dx).* (f_w(2:Ny,1).* Vx(2:Ny,2)) -(1/dy).*(f_w(2:Ny,1).* Vy(3:Ny+1,1)-f_w(1:Ny-1,1).* Vy(2:Ny,1)); % E,N,S
 R(1,2:Nx) = q(1,2:Nx) - (phi./dt).*(Sw_nu(1,2:Nx) - Sw_n(1,2:Nx)) - (1/dx).* (f_w(1,2:Nx).* Vx(1,3:Nx+1)-(f_w(1,1:Nx-1).* Vx(1,2:Nx))) -(1/dy).*(f_w(1,2:Nx).* Vy(2,2:Nx)); % E,W,N
% Residual of remaining cells 
%% how to define this in 2D
for i = 2:Ny
    for j = 2:Nx 
        R(i,j) = q(1,i*j) - (phi/dt).*(Sw_nu(i,j) - Sw_n(i,j)) - (1/dx).* (f_w(i,j).* Vx(i,j+1) - f_w(i,j-1).*Vx(i,j)) -(1/dy).* (f_w(i,j).* Vy(i+1,j) - f_w(i-1,j).*Vy(i,j));
    % check the above step ????????
    end
end 

end 
