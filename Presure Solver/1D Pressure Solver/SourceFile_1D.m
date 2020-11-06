                    %% 1D PRESSURE SOLVER 
%In compressible, slightly compressible and Fully compressible


  
clear all; close all; clc;
%% Initialization, Parameters, Grid
L            = 1 ;                      %length of the reservoir [m]
N            = 20;                      %number of grid cells
dx           = L/N;                     %Grid size
x            = linspace(dx/2,L-dx/2,N); %Location of the cell centre
xi           = linspace(0,L,N+1);       %Location of interfaces
Lamda        = zeros(N,1);              %initialization of vector

%% 
H = 1;  
%H = 0 for heterogeneous reservoir
%H = 1 for homogeneous reservoir

switch H
    case 0                              % 0 for heterogeneous reservoir
      Lamda(1:N) = 100.*(10.^rand(N,1));
    case 1                              % 1 for homogeneous reservoir
     Lamda(1:N)  = 1;   
end

Lamda=Lamda';                           % Transposing the matrix
%%
PL           = 1;                       %boundary condition at left boundary
PR           = 0;                       %boundary condition at right boundary
ceff=.01;
phi=.2;
p=zeros(N,1);
dt=0.00001;

%% Parameters needed for  compressible flow

rhoo = 1;           % Initial Density
phio = 0.3;         % Initial Porosity
po = 0;             % Initial Pressure
c_rho = 1;          % Initial Compressibility
c_phi = 1;
cell = [1, N];

pw = [2 0];    
PI = [1000 1000];
%%
[T, LamdaH] = Transmissibility(N,dx,Lamda); % Transmissibility Function
%%
flow = 3;
% flow = 1, Incompressible
% flow = 2, Slightly_compressible
% flow = 3, Fully compressible

switch flow 
    %% Incompressible 
    case 1
        % A*p = q
        [T, LamdaH] = Transmissibility(N,dx,Lamda);
                p = zeros(N,1);
                q = zeros(N,1);
%% ANALYTICAL SOLUTION

Pana = (PR - PL)*(x/L) + PL;
%%
% NUMERICAL SOLUTION
               [A] = compute_A(N, T); %calling function which computes A

        % Adding Wells
                wells = 1;          %0 = No wells, 1 = Wells
                
        switch wells 
             case(0)                % No wells, 
                  i= 1;
                      A(i,i) = A(i,i)+ T(i+1);
                      q(i) = q(i)*PL;
                  i = N;
                      A(i,i) = A(i,i)+ T(i+1);
                      q(i) = q(i)*PR;

             case(1)
                      pw = [1 0];   % No flow, wells are present
                      PI = [1000 1000];
                      cell = [1,N];
  
            [A,q] = addwells(A,q,Lamda,pw,PI,cell);
        end  

        figure (1)
        hold on 
        p = A\q;
        plot (x,p);
        xlabel('Length'); ylabel('Pressure');

        plot(x, Pana)           % Plotting analytical solution
        hold off
    title ('Numerical & Analytical Pressure solver');
    legend(' Analytical solution', 'Numerical solution');
%% Generating Error plot
Pana = Pana';
 E = (p - Pana);
 figure (3)
 plot (x, E,'Marker','*')
 xlabel('x'); ylabel('Error');
 title ('Error,1D Incompressible');

%% Slightly Compressible 
    case 2
      % best graph can be seen for dt = 0.00001 
    for ndt=1:100
    
    q=zeros(N,1);
    C=(phi*ceff/dt)*eye(N);
    PI = [1000,1000];
    pw=[1,0];
    cell=[1,N];
    
  [T, LamdaH] = Transmissibility(N,dx,Lamda);
  
    [A] = compute_A(N, T);
    [A,q] = addwells(A,q,Lamda,pw,PI,cell);
    
    implicit = 1;
    switch implicit
        case 0                              %   Explicit Solution
            p = C\(q+(C*p-A*p));           
        case 1                              %   Implicit solution
            p=(C+A)\(q+C*p);               
    end 
    
    [V] = ComputeVelocity(LamdaH,dx,N, p, PL, PR);
    v(:,ndt) = V;
    
    plot(x, p);
    hold on;
    
    end

    title ('Pressure Solver for Sli Compressible')
    xlabel('space'); ylabel('Pressure');
    
    figure (2)   % Velocity Plot
    plot (xi ,v);
    hold on;
    xlabel ('xi'); ylabel('velocity');
    title ('Velocity solver for Sli Compressible')
    
    compute_pana();

    
   %%  Fully compressible 
    case 3
        % best graph is seen for dt = 0.01 
         p = zeros(N,1);
         q = zeros(N,1);
        
        % Calling porosity, density and Transmisibility functions
        [phi, dphidp] = ComputePorosity(p, po,phio,c_rho);
        [rho, drhodp] = ComputeDensity(p, po,rhoo,c_phi);
        [T,LamdaH]    = Transmissibility(N,dx,Lamda.*rho);
        
        q = computewellsfluxes(pw, PI, Lamda, cell, p);

for ndt=1:200
    %Newton loop
    converged = 0;   % What is the significance of this ?
    rho_n = rho;     % are these density at nth time step ?
    phi_n = phi;
    %find p^n+2
    Residual = ComputeResidual(rho_n, phi_n, rho, phi, T, dt, p, q, N);  %Compute Transm
     iter_step = 1;
    while converged == 0 && iter_step < 100
        %Build Jacobian matrix
        
        [A] = compute_A(N, T); %calling function which computes A

        % V/dt(phi*rho)
        %C(i,i) = V/dt(dphidp*rho(i) +
        vec = (1/dt)*(dphidp .* rho + drhodp .* phi);
        C = diag(vec);
        %% Wells
        %q = rho*Lamda*PI*(Pwell - p_i) pi = pcell (RATES produced or injected)
        %W(i,i) = Lamda*PI(i)*rho(i);
        
        W =  zeros(N,N);  % this is diagonal matrix with the above in diagonal elements
        for i = 1:length(PI)         % computing w for all the wells
            W(cell(i), cell(i)) =  W(cell(i), cell(i))+ Lamda(cell(i))*PI(i)*rho(cell(i)) ;    % location at cell i  denisty at location
        end
        
        J = A + C + W;              % what is W ?????????
        dp = J\Residual;
        % Find pressure nu+1
        p = p+dp;     % increasing the pressure by dp
        
        % updating rock and fluid parameters
        [phi, dphidp] = ComputePorosity(p, po,phio,c_rho);
        [rho, drhodp] = ComputeDensity(p, po,rhoo,c_phi);
        T = Transmissibility(N , dx, Lamda.*rho);
        
        q = computewellsfluxes(pw, PI, Lamda, cell, p);
        
        % compute Residual
        Residual = ComputeResidual(rho_n, phi_n, rho, phi, T, dt, p, q, N); % This variables must match with the function variables
        iter_step  = iter_step + 1;
        
        if norm(Residual, inf) < 1e-6   % in place of infinty, you give 2nd norm which is the average square root of error at all the locations
            converged = 1;
        end
    end
        [V] = ComputeVelocity(LamdaH,dx,N, p, PL, PR); % Velocity
         v(:,ndt) = V;
         
        plot (x, p);
        hold on;
end
    xlabel('x');
    ylabel('Pressure');
    title('Pressure Solver');
    
    figure (2)
    plot (xi ,v);
    hold on;
    xlabel ('xi'); ylabel('velocity');
    title ('Velocity solver for Compressible')
    
end 
  