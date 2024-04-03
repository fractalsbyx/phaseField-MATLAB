%% periodic boundary condition
function dArrdx = gradArr_1D(Arr, rho, dx)
  drhodx_central  = (circshift(rho, -1) - circshift(rho, 1)) / 2.0 / dx;

  dArrdx_forward  = (circshift(Arr, -1) - Arr) / dx;
  dArrdx_backward = (Arr - circshift(Arr, 1)) / dx;
  
  dArrdx = zeros(size(drhodx_central));

  dArrdx(drhodx_central >= 0) = dArrdx_backward(drhodx_central >= 0);
  dArrdx(drhodx_central < 0) = dArrdx_forward(drhodx_central < 0);
end