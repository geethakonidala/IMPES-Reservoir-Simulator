                        %% Saturation Transport Solver
 %                           Geetha Krishna-4627628 
 
%% Retrieving velocity from Flow solver

clear all;close all; clc;

L       = 10;
N       = 200;                   % Number of Grid Cells

PL = 1;                          % Input BC
PR = 0;                          % Output BC
dx = L/N;                        % Grid Size
x = linspace(dx/2, L-dx/2,N);    % Location of Grid Center 
xi = linspace(0,L,N+1);          % Location of Interface

Lamda   = zeros(N,1);            % Initialization Vector

Homogeneous = 1;                
switch Homogeneous
    case 0                              % 0 for heterogeneous reservoir
      Lamda(1:N) = 100.*(10.^rand(N,1));
    case 1                              % 1 for homogeneous reservoir
     Lamda(1:N)  = 1;   
end

Lamda = Lamda';

[T,LamdaH] = Transmissibility(N,dx,Lamda);
 %% NUMERICAL PRESSURE SOLVER
 
    p = zeros(N,1);
    q = zeros(N,1);

[A] = compute_A(N, T);
%% Adding Wells

wells = 0;
switch wells 
    case(0)
        i= 1;
            A(i,i) = A(i,i)+ T(i+1);
            q(i) = q(i)*PL;
        i = N;
            A(i,i) = A(i,i)+ T(i+1);
            q(i) = q(i)*PR;

 case(1)
             pw = [1 0];
             PI = [1000 1000];
             cell = [1,N];
        %[A,q] = addwells(A,q,Lamda,pw,PI,cell)
            [A,q] = addwells(A,q,Lamda,pw,PI,cell);
end

p = A\q;
%% Initailizing parameters for Transport Solver
phi     = 0.3;      
nw      = 3.5;
no      = 2;
swc     = 0.2;
sor     = 0.1;
krwe    = 0.7;
kroe    = 0.8;
phi     = 0.3;
mu_w    = 10^-3;
mu_o    = 10^-2; 

Lamda_W = zeros(1, N);
Lamda_O = zeros(1, N);
q       = zeros(1, N);
Sw      = swc.*ones(1, N);
f_w     = zeros (1, N);

% Corey Relative Permeability Relationship
Krw     = krwe*((Sw-swc)/(1-swc-sor)).^nw;
Kro     = kroe*((1-Sw-sor)/(1-swc-sor)).^no;
                          %% Velocity Vector
Retrieved = 0;
% Retrieved = 1; Velocity extracted from Flow Solver
% Retrieved = 0; Velcocity assigned manually 
    switch Retrieved
        case 1
             v = zeros(1,N+1);

            v(1) = -Lamda(1)*2*(p(1)-PL)/dx;
            v(N+1) = -Lamda(N)*2*(PR-p(N))/dx;

            for i = 2:N
                v(i) = -LamdaH(i)*(p(i)-p(i-1))/dx;
            end 

         case 0
            v       = ones(1,N+1);
            v(1)    = 0; v(N+1) = 0;
    end 
    %% Flux 
     q(1)   = 1; q(N) = -1;
     f_inj  = 1;            % Water is injected at the left interface
%% 
Solver = 0;
% Solver = 0 for Explicit Solver 
% Solver = 1 for Implicit Solver
%%                              EXPLICIT SOLVER
switch Solver 
    case 0                     
 
    [f_w, dfw_ds ] = Compute_fractnlflow(mu_w, mu_o, Sw,krwe,kroe,sor,swc, nw, no, N);

    Sw(1) = 1;      % Saturation at injector
    dt= 0.001;
    ndt = 500;
    t = dt*ndt;

      %_______________________________________________________________%
 
for ndt = 1:ndt

           [f_w, dfw_ds ] = Compute_fractnlflow(mu_w, mu_o, Sw,krwe,kroe,sor,swc, nw, no, N);
       for j = 2:N-1                               
            Sn(j) = q(j)*dt/phi + Sw(j) - (dt ./(dx.* phi)).*(f_w(j).* v(j+1)- f_w(j-1).* v(j));
            Sw(j) = Sn(j);
       end 
            figure(1)
            hold on
            plot (x, Sw) 
end
    xlabel('x[mtr]'); ylabel(' Water Saturation, Sw');
    title (' Saturation Profile');

%% Analytical and Numerical Solution Comparison

        [x_ana, sw2] = Compute_Ana(t,no,nw,mu_w,mu_o,sor,swc,krwe,kroe,phi);   % Calling Analytical Solution Function
figure(2)
        hold on
        plot (x, Sw)
        plot(x_ana,sw2,'.-','linewidth',2)
        legend('Numerical Solution', 'Analytical Solution');
        xlabel('x[mtr]');  ylabel(' Water Saturation, Sw');
        title('Analytical Vs Numerical Solution');
 
%%Error Analysis
        sw2_i   = sw2(3:end);
        x_ana_i = x_ana(3:end);
        r       = round(N./3);
        x_b     = x(1:r);
        Sw_i    = interp1(x_ana_i,sw2_i,x_b,'spline');
        Sw_a = Sw(1:r);
        Error = max(abs(Sw_a - Sw_i));   % Maximum Error for that Resolution
 
        %%                  IMPLICIT SOLVER
    case 1                   

        Sw      = swc.* ones(N,1);
        Sw(1)   =1-sor;
        Krw     = zeros(N,1);
        Kro     = zeros(N,1);
        V       = zeros(N,N);
        q       = zeros(1,N);
        q(1)    = 1/dx; 
        q(N)    = -1/dx;

% relperms + mobility + fracflow
        Krw   = krwe*((Sw-swc)/(1-swc-sor)).^nw;
        Kro   = kroe*((1-Sw-sor)/(1-swc-sor)).^no;
%% 
        dt= 0.001;   % Time Step
        ndt = 500;   % Number of Time Steps
        t = dt*ndt;
        Sw(1)   =1-sor;
for j = 1:ndt
    
    Converged = 0;
    Sw_new = Sw; % saturation at nth step
    [f_w, dfw_ds ] = Compute_fractnlflow(mu_w, mu_o, Sw,krwe,kroe,sor,swc, nw, no, N);
    
    
    %dfwds = 2*(Sw.*(1-Sw)./(mu_o .* mu_w))./ ((Sw.^2 ./mu_w) + (((1-Sw).^2) ./mu_o));
    Residual = Compute_Residual(q, phi, dx,dt, f_w, Sw_new,Sw, v, N);
   
    iter_step = 1;
    while Converged ==0 && iter_step < 100
    %% Constructing Jacobian  J =  V 
    V = zeros(N,N);
    %Computing viscous term 
    for j=1:N
        if (j>1)
            V(j, j - 1) =  -dfw_ds(j -1).* v(j)/dx;
        end 
        if j< N
            V(j, j) = phi./dt +  dfw_ds (j).* v(j + 1)/dx;

        end 
    end 
    
    V(N, N) = phi./dt +  dfw_ds (N).* v(N + 1)/dx;

% V(1,1) = phi./dt;
%     for j = 2:N
%          V(j,j-1) = -dfw_ds(j -1).* v(j)/dx;
%          V(j,j)   = phi./dt +  dfw_ds (j).* v(j + 1)/dx;
%     end
     % V(n,n) has only phi/dt, since velocity at the boundaries are zero
    J =   V - diag(min(q',0).* dfw_ds);
       
    dSw = J\Residual;
    Sw = Sw + dSw;   % New Saturation
    Sw(N,1) = swc;
    %UPDATING Saturation 
    [f_w, dfw_ds ] = Compute_fractnlflow(mu_w, mu_o, Sw,krwe,kroe,sor,swc, nw, no, N);
    
    Residual = Compute_Residual(q, phi, dx,dt, f_w, Sw_new,Sw, v, N);
    
    norm(Residual, inf);
    
        if norm(Residual, inf) < 10^-6   % in place of infinty, you give 2nd norm which is the average square root of error at all the locations
           converged = 1;
        end
        
      iter_step = iter_step +1;
      
    end
    hold on
    plot (x, Sw);
end 

 xlabel('x[mtr]','linewidth',1.5); ylabel('Sw', 'linewidth',1.5);
 title('Saturation Profile,Implicit solver','linewidth',2);
 
end 
