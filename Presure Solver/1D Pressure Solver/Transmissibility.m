function [T, LamdaH] = Transmissibility(N,dx, Lamda)
LamdaH = zeros(N+1,1);
T = zeros(N+1,1);
LamdaH(1) = Lamda(1);
LamdaH(N+1) = Lamda(N);
LamdaH(2:N) = 2*Lamda(2:N).*Lamda(1:N-1)./(Lamda(2:N) + Lamda(1:N-1));

T(1) = LamdaH(1)./(dx^2/2);
T(N+1) = LamdaH(N+1)./(dx^2/2);
T(2:N) = LamdaH(2:N)./(dx^2);
end 