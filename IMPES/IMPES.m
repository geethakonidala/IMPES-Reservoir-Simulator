%%                              IMPES
%%                   @GEETHA KRISHNA CHOUDARY- 4627628
%                             12-6-2017
%%
clear all; close all; clc;
% Input parameters
Lx      = 10;                         % x-dimension [m]
Ly      = 10;                         % y-dimension [m]
Nx      = 20;                         % Number of grids in X-axis
Ny      = 20;                         % Number of grids in y-axis
N       = Nx * Ny;                    % Total number of grids in 2D
dx      = Lx / Nx;                    % Grid Size in X-axis
dy      = Ly / Ny;                    % % Grid Size in X-axis
x       = linspace(dx/2,Lx-dx/2,Nx);  %Location of cell centre
xi      = linspace(0,Lx,Nx+1);        %Location of interfaces
y       = linspace(dy/2,Ly-dy/2,Ny); 
yi      = linspace(0,Ly,Ny+1); 
N       = Nx*Ny;

%%  Permeability distribution 
H = 1;  
%H = 1 for Homogeneous Reservoir
%H = 0 for Heterogeneous Reservoir
switch H
    case 0                              % 0 for heterogeneous reservoir
        K= (1e-12).*(10.^rand(Ny,Nx));
    case 1                              % 1 for homogeneous reservoir
        K= (1e-12).* ones(Ny, Nx);      % Permeability [m^2];   
end
%% Initializing for wells
pw      = [10^6 0];                     % well pressure     
PI      = [1000 1000];                  % Productivity Index
cell    = [1, Nx*Ny];                   % Well Location 

%% Initailizing parameters for Transport Solver

phi     = 0.3;                        % Porosity
nw      = 3.5;                        % Corey exponent for water
no      = 2;                          % Corey exponent for oil
swc     = 0.2;                        % Connate water saturation
sor     = 0.1;                        % Residual Oil Saturation
krwe    = 0.7;                        % End point Kr of Water
kroe    = 0.8;                        % End point Kr of Oil
mu_w    = 10^-3;                      % Water Viscosity
mu_o    = 10^-2;                      % Oil Viscosity

%% Fixing the saturation at injection
Sw      = swc.*ones(Ny, Nx);          %initial water saturation in reservoir.
Sw(1,1) = 1-sor;                      % Injector at 1st cell  always water entering in 

CFL = 0.7;
dt= 100;
Ndt = 100;
for ndt = 1:Ndt
    A       = zeros(Nx*Ny , Nx*Ny);               
    p       = zeros(Nx*Ny, 1);
    q       = zeros(Nx*Ny, 1);
    %% Relative permeability and mobility  
    Krw         = krwe*((Sw-swc)/(1-swc-sor)).^nw;
    Kro         = kroe*((1-Sw-sor)/(1-swc-sor)).^no;
    Lamda_w     = zeros(Ny,Nx);
    Lamda_w     = K.*Krw./mu_w; 
    Lamda_o     = zeros(Ny,Nx);
    Lamda_o     = K.*Kro./mu_o;
    dLamda_w    =(K./mu_w)*(1/(1-swc-sor))*nw*krwe.*((Sw-swc)./(1-swc-sor)).^(nw-1);
    dLamda_o    =(K./mu_o)*(-1/(1-swc-sor))*no*kroe.*((1-Sw-sor)./(1-swc-sor)).^(no-1);
    lamda_total = Lamda_w+Lamda_o;
    
    %% Updating transimissibility 
    [Tx,Ty, LamdaHx, LamdaHy] = Trans2D_IMPES(dx,dy,lamda_total,Nx,Ny);
    
    %% Pressure solver
    [A] = computeA_2D(Nx, Ny, Tx, Ty);
    [A,q] = addwells(A,q,lamda_total,pw,PI,cell);
    p = A\q;
    P = reshape(p, Nx, Ny);
    [Vx,Vy] = ComputeVelocity_2D(LamdaHx, LamdaHy, Nx, Ny, P);
    [q] = ComputeWellsFluxes(P,PI,lamda_total,cell,pw,N); 
    q=reshape(q,Ny,Nx);  
    
    %% dt using CFL Condition
    if ndt ==1
        dfw_ds  = Computedfw_ds(Nx, Ny); 
       [dt] = Compute_CFL(phi,dx,dfw_ds,CFL,Vx,Vy, Nx, Ny);
    end 
     
    %% Solving for Saturation
    % fractioanl flow
    [f_w, dfw_ds] = Compute_fracflow_2D(Lamda_w,Lamda_o,dLamda_w,dLamda_o); 
    for i=1:Nx
        for j=1:Ny
   %Left boundaries
            if i==1&&j>1  
                sw_n(i,j)=dt*q(i,j)/phi-(dt/phi/dx)*(f_w(i,j)*Vx(i+1,j))-(dt/phi/dy)*(f_w(i,j)*Vy(i,j+1)-f_w(i,j-1)*Vy(i,j))+Sw(i,j);
            end
   %Bottom boundaries
            if j==1&&i>1  
                sw_n(i,j)=dt*q(i,j)/phi-(dt/phi/dx)*(f_w(i,j)*Vx(i+1,j)-f_w(i-1,j).*Vx(i,j))-(dt/phi/dy)*(f_w(i,j)*Vy(i,j+1))+Sw(i,j);  
            end
            if i>1&&j>1
   % Cells excluding boundaries
                sw_n(i,j)=dt*q(i,j)/phi-(dt/phi/dx)*(f_w(i,j)*Vx(i+1,j)-f_w(i-1,j).*Vx(i,j))-(dt/phi/dy)*(f_w(i,j)*Vy(i,j+1)-f_w(i,j-1)*Vy(i,j))+Sw(i,j);
            end
        end
    end
    
    sw_n(1,1) = 1-sor; 
    sw_n(Ny,Nx) = swc;
    figure (2)
    Sw = sw_n; 

end
figure (1)
    surf (P)
    xlabel(' Length , x'); 
    ylabel(' Length , y'); 
    zlabel(' Water Saturation, Sw');    
    title (' IMPES Pressure profile ');
 figure (2)
   surf (Sw)
   xlabel(' Length , x'); 
   ylabel(' Length , y'); 
   zlabel(' Water Saturation, Sw');
   title (' IMPES Saturation Profile  ');
