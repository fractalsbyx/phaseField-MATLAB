function D = div(vx, vy, rho, dx, dy)
  [xx, ~] = gradArr_2D(vx, rho, dx, dy);
  [~, yy] = gradArr_2D(vy, rho, dx, dy);
  D = xx+yy;
end