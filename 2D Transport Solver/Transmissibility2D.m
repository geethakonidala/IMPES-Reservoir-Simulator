function [Tx,Ty, LamdaHx, LamdaHy] = Transmissibility2D(Nx,Ny, dx,dy, Lamda)
% Computing Lamda & Transmissibility
%Initializing
LamdaHx = zeros(Nx +1,Ny);
LamdaHy = zeros(Nx, Ny+1);
Tx      = zeros(Nx+1,Ny);
Ty      = zeros(Nx, Ny +1);

%Harmonic Averaging 
LamdaHx(1,:)     = Lamda(1,:);    %  values at boundaries  
LamdaHx(2:Nx,:)  = 2*Lamda(2:Nx,:).*Lamda(1:Nx-1,:)./(Lamda(2:Nx,:) + Lamda(1:Nx-1,:));
LamdaHx(Nx+1,:)  = Lamda(Nx,:);

LamdaHy(:,1)     = Lamda(:,1);
LamdaHy(:,2:Ny)  = 2*Lamda(:,2:Ny).*Lamda(:,1:Ny-1)./(Lamda(:,2:Ny) + Lamda(:,1:Ny-1));
LamdaHy(:,Ny+1)  = Lamda(:,Ny);

%Transmissibility averaging 
Tx(1,:)          = LamdaHx(1,:)./(dx^2/2);
Tx(Nx+1,:)       = LamdaHx(Nx+1,:)./(dx^2/2);
Tx(2:Nx,:)       = LamdaHx(2:Nx,:)./(dx^2);

Ty(:,1)          = LamdaHy(:,1)./(dy^2/2);
Ty(:,Ny+1)       = LamdaHy(:,Ny+1)./(dy^2/2);
Ty(:,2:Ny)       = LamdaHy(:,2:Ny)./(dy^2);

end 