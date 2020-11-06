function q = computewellsfluxes2D(pw, PI, Lamda, cell, p, Nx, Ny)

    q = zeros(Nx*Ny,1);   %Need to define Nx & Ny 

for w = 1:length(PI)
    
    q(cell(w)) = q(cell(w),1)+PI(w)*Lamda(cell(w))*(pw(w) - p(cell(w)));   % check this line
end
end