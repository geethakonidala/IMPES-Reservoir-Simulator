function dt_s = Compute_Syability(CFL,dfw_ds,phi,dx)
dfw_dsm = max(dfw_ds);
        dt_s = CFL.*dx.*phi./(dfw_dsm);
end 