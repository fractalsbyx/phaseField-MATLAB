%% periodic boundary condition
function [dArrdx, dArrdy] = eno_gradArr_2D(Arr, rho, dx, dy)
  drhodx_central  = (circshift(rho, [0, -1]) - circshift(rho, [0, 1])) / 2.0 / dx;
  drhody_central  = (circshift(rho, [-1, 0]) - circshift(rho, [1, 0])) / 2.0 / dx;

  dArrdx_forward  = (circshift(Arr, [0, -1]) - Arr) / dx;
  dArrdx_backward = (Arr - circshift(Arr, [0, 1])) / dx;
  dArrdy_forward  = (circshift(Arr, [-1, 0]) - Arr) / dy;
  dArrdy_backward = (Arr - circshift(Arr, [1, 0])) / dy;
  
  dArrdx = zeros(size(drhodx_central));
  dArrdy = zeros(size(drhody_central));
  
  dArrdx(drhodx_central >= 0) = dArrdx_backward(drhodx_central >= 0);
  dArrdx(drhodx_central < 0) = dArrdx_forward(drhodx_central < 0);
  dArrdy(drhody_central >= 0) = dArrdy_backward(drhody_central >= 0);
  dArrdy(drhody_central < 0) = dArrdy_forward(drhody_central < 0);
end