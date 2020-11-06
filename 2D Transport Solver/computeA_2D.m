function [A] = computeA_2D(Nx, Ny, Tx, Ty)

A             = zeros(Nx*Ny , Nx*Ny);  

for i = 1:Nx
    for j =1:Ny
        
        % Indexing the cells around i
        I  = (j-1)*Nx  +  i;            % cell i 
        Ie = (j -1)*Nx + (i+1);         % east
        Iw = (j-1)*Nx  + (i-1);         % west 
        In = (j*Nx)    +  i;            % north
        Is = (j-2)*Nx  +  i;            % south
        
        if i > 1            % adding West connection
            A(I,I) = A(I,I) + Tx(i,j);         % Interface between west and i is i
            A(I,Iw) = -Tx(i,j);
        end 
        
        if i < Nx           % adding East connection 
            A(I,I) = A(I, I) +  Tx(i+1,j);     % Interface between i and east is i+1
            A(I,Ie)= -Tx(i+1,j);
        end 
        
        if j > 1            % adding South connection
           A(I,I)  = A(I, I)  + Ty(i,j);
           A(I,Is) = -Ty(i,j);       % Interface between south and i is j
        end 
        
        if j < Ny         % adding North connection
            A(I,I) = A(I, I)  + Ty(i,j+1);
            A(I,In) = -Ty(i,j+1);     % Interface between i and north is j+1
        end 
    end
end 