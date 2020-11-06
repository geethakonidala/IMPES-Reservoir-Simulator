function [f_w, dfw_ds] = Compute_fracflow_2D(lamda_w,lamda_o,dlamda_w,dlamda_o)
  f_w    =lamda_w./(lamda_w+lamda_o);
  dfw_ds = (dlamda_w.*(lamda_w+lamda_o)-lamda_w.*(dlamda_w+dlamda_o))./(lamda_w+lamda_o).^2;
end 