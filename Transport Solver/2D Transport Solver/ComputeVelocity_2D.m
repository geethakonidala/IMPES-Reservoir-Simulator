function [ Vx, Vy] = ComputeVelocity_2D(LamdaHx, LamdaHy, Nx, Ny, P, dx, dy)

LamdaHx = LamdaHx';
LamdaHy = LamdaHy';
% P = P';
%Calculating velocity, exculding the boundaries
Vx = zeros(Ny, Nx+1);  % storing x in columns, y in rows
Vy = zeros(Ny+1,Nx);

%Velocity in X Direction

for i = 1:Ny
    Vx(i,2:Nx) = -LamdaHx(i,2:Nx).*(P(i,2:Nx)- P(i,1:(Nx-1)))/dx ;
end 

%Velocity in Y Direction

for i = 1:Nx
    Vy(2:Ny,i) = -LamdaHx(2:Ny,i).*(P(2:Ny,i)- P(1:(Nx-1),i))/dy ;
end 
end 

