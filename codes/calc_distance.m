function [dist] = calc_distance(x1, y1, x2, y2, lx, ly)
  dist = sqrt(periodic_diff(x1,x2,lx).^2 + periodic_diff(y1,y2,ly).^2);
end