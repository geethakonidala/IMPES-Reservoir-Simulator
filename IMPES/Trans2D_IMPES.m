function [Tx,Ty, LamdaHx, LamdaHy] = Trans2D_IMPES(dx,dy,lamda_total,Nx,Ny)

Tx      = zeros(Ny,Nx+1);       % Indexing x axis in columns
Ty      = zeros(Ny +1, Nx);     % y axis in rows

LamdaHx=zeros(Nx+1,Ny);        %Harmonic avg in x
LamdaHy=zeros(Nx,Ny+1);        %Harmonic avg in y

LamdaHx(1,1:Ny)=lamda_total(1,1:Ny);        %left bound
Tx(1,1:Ny)=LamdaHx(1,1:Ny)./(dx^2/2);

LamdaHx(Nx+1,1:Ny)=lamda_total(Nx,1:Ny);    %right bound
Tx(Nx+1,1:Ny)=LamdaHx(Nx+1,1:Ny)./(dx^2/2);

LamdaHy(1:Nx,1)=lamda_total(1:Nx,1);        %bottom bound
Ty(1:Nx,1)=LamdaHy(1:Nx,1)./(dy^2/2);
 
LamdaHy(1:Nx,Ny+1)=lamda_total(1:Nx,Ny);    %top bound
Ty(1:Nx,Ny+1)=LamdaHy(1:Nx,Ny+1)./(dy^2/2);

%% Tx 
% i for Nx j for Ny
for i=2:Nx
    for j=1:Ny
        LamdaHx(i,j)=2*lamda_total(i-1,j)*lamda_total(i,j)/(lamda_total(i-1,j)+lamda_total(i,j));
        Tx(i,j)=LamdaHx(i,j)/(dx^2);
    end
end
%% Ty 
for i=1:Nx
    for j=2:Ny
        LamdaHy(i,j)=2*lamda_total(i,j-1)*lamda_total(i,j)/(lamda_total(i,j-1)+lamda_total(i,j));
        Ty(i,j)=LamdaHy(i,j)/(dy^2);
    end
end

end 
