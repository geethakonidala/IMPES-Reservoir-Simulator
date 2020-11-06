function [ Vx,Vy] = ComputeVelocity_2D(LamdaHx, LamdaHy, Nx, Ny, P)
%% Velocity in X Direction
Vx(1,1:Ny)=0;       %left
Vx(Nx+1,1:Ny)=0;    %right
for i=2:Nx
    for j=1:Ny
        Vx(i,j)=-LamdaHx(i,j)*(P(i,j)-P(i-1,j));
    end
end

%% Velocity in Y Direction
Vy(1:Nx,1)=0;     %bot
Vy(1:Nx,Ny+1)=0;     %top
for i=1:Nx
    for j=2:Ny
        Vy(i,j)=-LamdaHy(i,j)*(P(i,j)-P(i,j-1));
    end
end

end 

