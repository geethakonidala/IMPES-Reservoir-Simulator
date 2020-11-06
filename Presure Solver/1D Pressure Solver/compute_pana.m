           %% ANALYTICAL SOLUTION FOR SLIGHTLY COMPRESSIBLE

%Function for validating Error term for Slightly Compressible

function  PA= compute_pana()

Lg            = 100 ;                            %length of the reservoir [m]
N            = 21;                             %number of grid cells
dx           = 2*Lg/N;                           %Grid size
x            = linspace(-Lg ,Lg , N); %Location of the cell centre
xi           = linspace(0,Lg,N+1);               %Location of interfaces
Lamda       = zeros(N,1);                       %initialization of vector
xa          = linspace(-Lg, Lg,N);

Homogeneous = 1;                
switch Homogeneous
    case 0                              % 0 for heterogeneous reservoir
      Lamda(1:N) = 100.*(10.^rand(N,1));
    case 1                              % 1 for homogeneous reservoir
     Lamda(1:N)  = 1;   
end

Lamda=Lamda';


PL           = 0;                       %boundary condition at left boundary
PR           = 0;                       %boundary condition at right boundary
ceff= 0.01;
phi=.3;
p=zeros(N,1);
p((N-1)/2 +1,1) = 15;
dt=0.00001;
[T, LamdaH] = Transmissibility(N,dx,Lamda); % Transmissibility Function


%_______________________________________________________________%
 
 
for ndt=1:200
    
    q=zeros(N,1);
    C=(phi*ceff/dt)*eye(N);
    PI = [1000,1000];
    pw=[0,0];
    cell=[1,N];
  
    [A] = compute_A(N, T);  % calling computeA function to calculate A
    [A,q] = addwells(A,q,Lamda,pw,PI,cell);
    
    implicit = 1;
    switch implicit
        case 0                              %   Explicit Solution
            p = C\(q+(C*p-A*p));           
        case 1                              %   Implicit solution
            p=(C+A)\(q+C*p);               
    end 
    
    
end


%% 
K = 10^-14;
mu = 10^-03;
t = 10;
Ceff = 10^-06;

xa = linspace(-Lg, Lg,N);

D = K/(mu*phi*Ceff);
S  = (2*D*t)^0.5;

pana = exp(-0.5*(xa.^2)/S.^2)/(S*sqrt(2*3.14));


figure (3)
hold on
 plot (x, p)
plot (x, pana)
hold off
legend('numerical solution','analytical solution');
xlabel('space'); ylabel('Pressure');

end 