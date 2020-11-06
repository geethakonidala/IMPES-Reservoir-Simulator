function [B] = computeB_2D(Nx, Ny, dx, dy, Vx, Vy, dfw_ds)

B      = zeros(Nx*Ny , Nx*Ny);  
%          V(j,j-1) = -dfw_ds(j -1).* u(j)/dx;
%          V(j,j)   =   dfw_ds (j).* u(j + 1)/dx;

%         I  = (i-1)*Nx  +  j;            % cell i 
%         Ie = (i -1)*Nx + (j+1);         % east
%         Iw = (i-1)*Nx  + (j-1);         % west 
%         In = (i*Nx)    +  j;            % north
%         Is = (i-2)*Nx  +  j;            % south
        
        B(1,1)     = (1/dx).* dfw_ds(1,1).* Vx(1,2) + (1/dy).*dfw_ds(1,1).* Vy(2,1); %N,E
        B(2:Ny,1)  = (1/dx).* (dfw_ds(2:Ny,1).* Vx(2:Ny,2)) + (1/dy).*(dfw_ds(2:Ny,1).* Vy(3:Ny+1,1) - dfw_ds(1:Ny-1,1).* Vy(2:Ny,1)); %E,N,S
        B(1,2:Nx)  = (1/dx).* ((dfw_ds(1,2:Nx).* Vx(1,3:Nx+1))- (dfw_ds(1,1:Nx-1).* Vx(1,2:Nx))) + (1/dy).*dfw_ds(1,2:Nx).* Vy(2,2:Nx); %E,W,N
             %error vostey pai line lo check chey
 for j = 2:Nx
    for  i=2:Ny
         B(i,j) = (1/dx).* ((dfw_ds(i,j).* Vx(i,j+1))- (dfw_ds(i,j-1).* Vx(i,j))) + (1/dy).*(dfw_ds(i,j).* Vy(i+1,j)-dfw_ds(i-1,j).* Vy(i,j)); % E,W,N,S
%          B(Ny,2:Nx-1) =(1/dx).* ((dfw_ds(Ny,2:NX-1).* Vx(Ny,3:Nx))- (dfw_ds(Ny,1:Nx-2).* Vx(1,2:Nx-1))) + (1/dy).*(dfw_ds(Ny,2:Nx-1).* Vy(Ny+1,2:Nx-1)-dfw_ds(Ny-1,2:Nx-1).* Vy(Ny,2:Nx-1)); % E,W,N,S   
        
%         
%         if j > 1            % adding West connection
%             B(I,Iw) = B(I,I) - dfw_ds(i,j -1).* Vx(i,j)/dx ;    % Interface between west and i is i
% 
%         end 
%         
%         if j < Nx           % adding East connection 
%             B(I,Ie) = B(I, I) - dfw_ds(i,j).* Vx(i,j+1)/dx ;     % Interface between i and east is i+1
%         end 
%         
%         if i > 1            % adding South connection
%            B(I,Is) = B(I, I) - dfw_ds(i-1,j).* Vy(i,j)/dy;       % Interface between south and i is j
%         end 
%         
%         if i < Ny         % adding North connection
%             B(I,In) = B(I, I) - dfw_ds(i,j).* Vy(i+1,j)/dy;     % Interface between i and north is j+1
%         end 
    end
end 