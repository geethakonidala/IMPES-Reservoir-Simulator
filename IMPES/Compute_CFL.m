function [dt] = Compute_CFL(phi,dx,dfw_ds,CFL,Vx,Vy, Nx, Ny)
      
  Gx = Vx(2:Ny+1,:).*dfw_ds(:,:)./phi;
  Gy = Vy(:,2:Nx+1).*dfw_ds(:,:)./phi;
  
  a = max(Gx);
  c = max(a);
  
  b = max(Gy);
  d = max(b);
  
  dt = CFL*dx./max(c, d);
 
end 
