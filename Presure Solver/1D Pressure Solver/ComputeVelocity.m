function [V] = ComputeVelocity(LamdaH,dx,N, p, PL, PR)

V = zeros(N+1,1);

% for no flow condition
    for i = 2:N
        V(i)   = -LamdaH(i)* (p(i)-p(i-1))/dx;
       % V(1)   = -LamdaH(1)* 2*(p(1)-PL)/dx;
        %V(N+1) = -LamdaH(N)* 2*(PR - p(N))/dx;
    end
    
end 