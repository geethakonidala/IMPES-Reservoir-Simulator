%%  
%%Retrieving Velocity from flow solver

clear all; close all; clc;
%Initializing Parameters
Lx          = 10 ;                      
Ly          = 10 ; 
Nx          = 10;                      
Ny          = 10;

dx          = Lx/Nx;                     
dy          = Ly/Ny;

x           = linspace(dx/2,Lx-dx/2,Nx); %Location of cell centre
xi          = linspace(0,Lx,Nx+1);       %Location of interfaces
y           = linspace(dy/2,Ly-dy/2,Ny); 
yi          = linspace(0,Ly,Ny+1);

 
Lamda        = zeros(Nx,Ny);            %Initialization of vector
Lamda(1:Nx,:)  = 1;                     %Input lambda for homogeneous sol
Lamda(:,1:Ny)  = 1;
% Lamda(Nx,Ny) = 100.*(10.^rand(Nx,1)); %For hetrogeneous reservoir

rhoo    = 1;           % Initial Density
phi_o   = 0.3;         % Initial Porosity
po      = 0;           % Initial Pressure 
ceff    = 0.01;        % Initial Compressibility
c_rho   = 1;
c_phi   = 1;            
dt      = 0.01;       % Time Step

A = zeros(Nx*Ny , Nx*Ny);               
p = zeros(Nx*Ny, 1);
q = zeros(Nx*Ny, 1);

pw = [2 0];                             % Initializing for wells
PI = [1000 1000];
cell = [1, Nx*Ny];     

[Tx,Ty, LamdaHx, LamdaHy] = Transmissibility2D(Nx,Ny, dx,dy, Lamda);

% 
        [A] = computeA_2D(Nx, Ny, Tx, Ty);
        [A,q] = addwells(A,q,Lamda,pw,PI,cell);

        p = A\q; 
        P = reshape(p, Nx, Ny);
%         surf(P);     

% Velocity calculation
     Pa = P';
    [Vx, Vy] = ComputeVelocity_2D(LamdaHx, LamdaHy, Nx, Ny, Pa, dx, dy);


%% 2D Explicit solver

% parameters
mu_w        = 1.0e-3;
mu_o        = 10.0e-3;
nw          = 3.5;
no          = 2;
swc         = 0.2;
sor         = 0.1;
krwe        = 0.7;
kroe        = 0.8;
phi         = 0.3;

Lamda_W     = zeros(Ny,Nx);         % saving x axis in columns
Lamda_O     = zeros(Ny,Nx);
q           = zeros(1,Nx*Ny);       % Matrix dimension match avvakuntey ikkada chudu
Sw          = swc.*ones(Ny, Nx); 
Vx          = ones(Ny, Nx+1);       % storing x in columns, y in rows
Vy          = ones(Ny+1,Nx);
f_w         = zeros(Ny, Nx);
q(1,1)      = 1; q(1,Nx*Ny) = -1;   % assigning fluzes to 1st and last cell

%% Fixing the saturation at injection
Sw(1,1) = 1-sor;     % Injector at 1st cell
Sw(Ny,Nx) = 1-swc;   % Producer at last cell
%% Velocity 

Retrieved = 1;
% Retrieved = 1; Velocity extracted from Flow Solver
% Retrieved = 0; Velcocity assigned manually 
    switch Retrieved
        case 1
             [Vx, Vy] = ComputeVelocity_2D(LamdaHx, LamdaHy, Nx, Ny, Pa, dx, dy);
         case 0
            Vx(:,1) = 0; Vx(:, Nx+1) = 0;
            Vy(1,:) = 0; Vy(Ny+1,:)  = 0;
    end 
%making velocities at interface zero
%%
% relperms + mobility + fracflow
Krw         = krwe*((Sw-swc)/(1-swc-sor)).^nw;
Kro         = kroe*((1-Sw-sor)/(1-swc-sor)).^no;


[f_w, dfw_ds] = Compute_fracflow_2D( Sw, swc,krwe,kroe, sor, mu_w, mu_o,nw,no, Nx, Ny);
 
Sw_nu      = Sw;
dt=0.001;
ndt = 100;
 
 for i=1:ndt
       converged     = 0;
       Sw_n          = Sw_nu ;      % Saturation of Previous step;
       [f_w, dfw_ds] = Compute_fracflow_2D( Sw, swc,krwe,kroe, sor, mu_w, mu_o,nw,no, Nx, Ny);
      
       R = Compute_Residual2D(q, phi, dx,dy,dt, f_w, Sw_n,Sw_nu, Vx,Vy, Nx,Ny);
       iter_step =1 ;
       % Newton loop to calculate pressure %
while converged ==0   
            
            %A = zeros(Nx*Ny, Nx*Ny);
            vec = (phi./dt);
            A = diag(vec);
          
           [B] = computeB_2D(Nx, Ny, dx, dy, Vx, Vy, dfw_ds);
       
        J =   B + A;
        dSw = J\Residual';
      
        Sw_nu = Sw_nu + dSw;
   
     %UPDATING Saturation 
   [f_w, dfw_ds] = Compute_fracflow_2D( Sw, swc,krwe,kroe, sor, mu_w, mu_o,nw,no, Nx, Ny);
    
    R = Compute_Residual2D(q, phi, dx,dy,dt, f_w, Sw_n,Sw_nu, Vx,Vy, Nx,Ny);
    
    norm(Residual, inf);
    
    if norm(Residual, inf) < 10^-6   % in place of infinty, you give 2nd norm which is the average square root of error at all the locations
%            norm(Residual, inf)
           converged = 1;
    end
       
    iter_step = iter_step+1;  
end
     Sw(:,i)  = Sw_nu;
     Sa = reshape(Sw, Nx, Ny);
     hold on
     surf(Sa);
 
 end

    zlabel(' Water Saturation, Sw');
    title (' Explicit Saturation Profile (20x20)');
 