function Residual = ComputeResidual(rho_n, phi_n, rho, phi, T, dt, p, q, N)
%rho_n at previous time step and rho at nu, we consider for cell i

% constrction of A Matrix
A = zeros(N,N);
 for i=1:N
    if (i>1) % there is a left interface
        %T(i)*(p(i) - p(i-1)
            A(i,i)  =A(i,i)+T(i);
            A(i,i-1)= -T(i);
    end
    if (i<N) % there is a right neighbour
        % T(i+1)*(p(i)-p(i+1) )
            A(i,i)  =  A(i,i)+ T(i+1);
            A(i,i+1)= -T(i+1);
    end
end 
    Residual = rho .* q - ((rho .* phi - rho_n .* phi_n)/dt) - A*p;

end 

% instead of doing this, you can write t(i)(p(i)-p(i-1))+ .... like you did
% in flow sheet 
