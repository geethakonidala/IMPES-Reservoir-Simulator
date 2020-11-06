function CFL = Compute_CFL(dx,phi,u,dfw_ds)
CFL = dx.*phi./(u.*dfw_ds);
end 