function L = lap_1D(n, dx)
L = (circshift(n, 1) + circshift(n, -1) - 2.d0 * n) / dx / dx;
end
