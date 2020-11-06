                        %2D PRESSURE SOLVER
              %Imcompressible, Slightly compressible & Compressible flow          
                          % GEETHA KRISHNA
      %%
 close all; clc;

%Initializing Parameters
Lx            = 1 ;                      
Ly            = 1 ; 
Nx            = 10;                      
Ny            = 10;

dx           = Lx/Nx;                     
dy           = Ly/Ny;

x            = linspace(dx/2,Lx-dx/2,Nx); %Location of cell centre
xi           = linspace(0,Lx,Nx+1);       %Location of interfaces
y            = linspace(dy/2,Ly-dy/2,Ny); 
yi           = linspace(0,Ly,Ny+1);       


Lamda        = zeros(Nx,Ny);            %Initialization of vector
Lamda(1:Nx,:)  = 1;                     %Input lambda for homogeneous sol
Lamda(:,1:Ny)  = 1;
%Lamda(Nx,Ny) = 100.*(10.^rand(Nx,1)); %For hetrogeneous reservoir

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

%% 
flow = 3;
% flow = 1, Incompressible Solution
% flow = 2, Slightly_compressible Solution
% flow = 3, Fully compressible Solution

switch flow 
    case 1

        [A] = computeA_2D(Nx, Ny, Tx, Ty);
        [A,q] = addwells(A,q,Lamda,pw,PI,cell);

        p = A\q; 
        P = reshape(p, Nx, Ny);
        surf(P);

        title('2D Incompressible Solver');
        xlabel('Length in x');ylabel('Length in y'); zlabel('Pressure');

%% Velocity calculation
     Pa = P';
[Vx, Vy] = ComputeVelocity_2D(LamdaHx, LamdaHy, Nx, Ny, Pa, dx, dy);

figure (2)
subplot(1,2,1);
surf(Vx);
title ('Vx, 2D Incom'); 
xlabel('Nx');ylabel('Ny'); zlabel('Vx'); 

subplot(1,2,2);
surf(Vy);
title ('Vy, 2D Incom'); 
xlabel('Nx');ylabel('Ny'); zlabel('Vy'); 


  %% for Slightly compressible  
    case 2
       % [Tx,Ty, LamdaHx, LamdaHy] = Transmissibility2D(Nx,Ny, dx,dy, Lamda);
 
        for ndt=1:200
    
             q      = zeros(Nx*Ny,1);
             C      = (phi_o*ceff/dt)*eye(Nx*Ny);
            [A]     = computeA_2D(Nx, Ny, Tx, Ty);
    
            [A,q]   = addwells(A,q,Lamda,pw,PI,cell);
    
            implicit = 1;
                switch implicit
                    case 0                              %   Explicit Solution
                         p = C\(q+(C*p-A*p));           
                    case 1                              %   Implicit solution
                         p=(C+A)\(q+C*p);               
                end 
    
                hold on;
                P = reshape(p, Nx, Ny);
                surf(P);
        end
        hold off 
% P = reshape(p, Nx, Ny);
% surf(P);
title('2D Slightly Compressible solver');
xlabel('Length in x');ylabel('Length in y'); zlabel('Pressure');


[Vx, Vy] = ComputeVelocity_2D(LamdaHx, LamdaHy, Nx, Ny, P, dx, dy);
figure (2)
subplot(1,2,1);
surf(Vx);
title ('Vx, 2D Sli com'); 
xlabel('Nx');ylabel('Ny'); zlabel('Vx'); 

subplot(1,2,2);
surf(Vy);
title ('Vy, 2D Sli com'); 
xlabel('Nx');ylabel('Ny'); zlabel('Vy'); 

%% for fully compressible
    case 3 
        
        q = computewellsfluxes2D(pw, PI, Lamda, cell, p,Nx,Ny);
        iter_step = 1;

        for ndt = 1:200
             %Newton loop
             converged = 0;  
             rho_n = rho;    
             phi_n = phi;
             %find p^n+2
             Residual = ComputeResidual2D(rho_n, phi_n, rho, phi, Tx,Ty, dt, p, q, Nx, Ny);  %Compute Transm
    
        while converged == 0 && iter_step < 100
   
             [A] = computeA_2D(Nx, Ny, Tx, Ty);
     
             vec = (1/dt)*(dphidp .* rho + drhodp .* phi);
             C   = diag(vec);

            %% Wells
        %q = rho*Lamda*PI*(Pwell - p_i) pi = pcell (RATES produced or injected)
        %W(i,i) = Lamda*PI(i)*rho(i);
        
             W =  zeros(Nx*Ny,Nx*Ny);  % diagonal matrix with the above in diagonal elements
        
        for i = 1:length(PI)      % Either PI/pw. this iteration is from 1 to nbr of wells
                W(cell(i), cell(i)) =  W(cell(i), cell(i))+ Lamda(cell(i))*PI(i)*rho(cell(i)) ;    % location at cell i  denisty at location
        end
        
            J = A + C + W;             
            dp = J\Residual;
        % Finding pressure pressure for next iteration
            p = p+dp;     % increasing the pressure by dp
        
        % update rock and fluid parameters
            [phi, dphidp] = ComputePorosity(p, po,phio,c_rho);
            [rho, drhodp] = ComputeDensity(p, po,rhoo,c_phi);
            rho_ma=reshape(rho,Nx,Ny);
            [Tx, Ty] =  Transmissibility2D(Nx,Ny, dx,dy, rho_ma*Lamda);
        
            q = computewellsfluxes2D(pw, PI, Lamda, cell, p, Nx, Ny);
        
        % updating Residual
            Residual = ComputeResidual2D(rho_n, phi_n, rho, phi, Tx,Ty, dt, p, q, Nx, Ny);  % This variables must match with the function variables
            max(Residual);
        
            Residualplot(iter_step) = max(Residual); % preallocation is not possible
            iter_step = iter_step +1;
        
             if norm(Residual, inf) < 1e-6   % in place of infinty, you give 2nd norm which is the average square root of error at all the locations
                    converged = 1;
             end
        end 
         
        figure (1)   
        hold on 
        P = reshape(p, Nx, Ny);
        surf(P);
       
        end
  hold off
  
xlabel('Length, x axis ');  ylabel('Length, y axis'); zlabel('Pressure')
title('2D Compressible solver');

[Vx, Vy] = ComputeVelocity_2D(LamdaHx, LamdaHy, Nx, Ny, P, dx, dy);
figure (2)
subplot(1,2,1);
surf(Vx);
title ('Vx, 2D Fully com'); 
xlabel('Nx');ylabel('Ny'); zlabel('Vx'); 

subplot(1,2,2);
surf(Vy);
title ('Vy, 2D Fully com'); 
xlabel('Nx');ylabel('Ny'); zlabel('Vy'); 

        
end 



