function [diff] = periodic_diff(x1, x2, lx)
  diff = (x1-x2);
  b = (abs(diff)>lx./2);
  diff = diff - sign(diff).*b.*(lx);
end