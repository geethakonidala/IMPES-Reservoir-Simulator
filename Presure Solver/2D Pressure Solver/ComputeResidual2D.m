function Residual = ComputeResidual2D(rho_n, phi_n, rho, phi, Tx,Ty, dt, p, q, Nx, Ny)

% Calling function which calculates A matrix
[A] = computeA_2D(Nx, Ny, Tx, Ty);  

Residual = rho .* q - ((rho .* phi - rho_n .* phi_n)/dt) - A*p;

end 