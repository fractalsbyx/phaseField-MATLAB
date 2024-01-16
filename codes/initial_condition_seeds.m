close all; clc; clear;
%% Setting up random seed
rng(2022);

%% Setting up simulation parameters
Params
dt = 0.1;

nprint = 100;
t0 = 0.d0;
thold = 9.0;

%% Generating grid points
[X, Y, Z] = meshgrid(x0:dx:xf, y0:dy:yf, z0:dz:zf);

%% Initializing order parameter n
n = zeros(Nx, Ny, Nz, np);
range = 4; % 3
coord = zeros(np, 3);

curr_id = 1;
exclude_distance = 25; %% 40 grid points
while curr_id <= np
  coord(curr_id, 1) = 1+floor(rand()*Nx);
  coord(curr_id, 2) = 1+floor(rand()*Ny);
  coord(curr_id, 3) = 1+floor(rand()*Nz);

  coord(1, 1) = 1+floor(Nx/2);
  coord(1, 2) = 1+floor(Ny/2);
  coord(1, 3) = 1+floor(Nz/2);
  
  is_valid_grain = true;
  for past_id = 1 : (curr_id - 1)
    distance = calc_distance(coord(past_id,:,:,:), coord(curr_id,:,:,:), Nx, Ny, Nz);
    if distance < exclude_distance
      fprintf("curr pos = [%i, %i], past pos = [%i, %i], distance = %f\n", coord(curr_id, 1), coord(curr_id, 2), coord(past_id, 1), coord(past_id, 2), distance);
      is_valid_grain = false;
      break;
    end
  end
  
  if is_valid_grain == true
    for x_id = (coord(curr_id, 1) - range) : (coord(curr_id, 1) + range)
      for y_id = (coord(curr_id, 2) - range) : (coord(curr_id, 2) + range)
        for z_id = (coord(curr_id, 3) - range) : (coord(curr_id, 3) + range)
          x_c  = 1+mod(x_id-1, Nx);
          y_c  = 1+mod(y_id-1, Nx);
          z_c  = 1+mod(z_id-1, Nx);
          if calc_distance([x_c, y_c, z_c], coord(curr_id,:,:,:), Nx, Ny, Nz) <= range
            n(x_c, y_c, z_c, curr_id) = 1.0;
          end
        end
      end
    end
    curr_id = curr_id + 1;
  end
end

%% Running simulations
t = t0;
tic;
while t < thold 
  toc;
  fprintf('Initializing, t = %f\n', t)
  
  % update eta (n)
  dndt = zeros(Nx, Ny, Nz, np);
  for istep = 1:nprint
    all_sum_nsq = sum(n.^2, 4);

    for i = 1:np
      sum_nsq = all_sum_nsq - n(:,:,:,i).^2; 
      term_sum_nsq = 2 * alpha * n(:,:,:,i) .* sum_nsq;
      % bulk term
      term1 = m0 * (-n(:,:,:,i) + n(:,:,:,i).^3 + term_sum_nsq);
      
      % interfacial term
      lap_n = lap(n(:,:,:,i), dx, dy, dz);
      term2 = -kappa * lap_n;

      % update dndt
      dndt(:,:,:,i) = -L * (term1 + term2);
    end
    
    % update n and rho
    n = n + dndt * dt;
    
    % update time t
    t = t + dt;
  end 
end



%% Write n to files
fileID = fopen('../data/initial_n.dat','w');
fwrite(fileID, n, 'double');
fclose(fileID);
fprintf("Initial condition written to ../data/initial_n.dat\n")

%% defining functions
function [diff] = periodic_diff(x1, x2, lx)
  diff = (x1-x2);
  b = (abs(diff)>lx/2);
  diff = diff - sign(diff).*b.*(lx);
end
function [dist] = calc_distanceB(x1, y1, x2, y2, z1, z2, lx, ly, lz)
  dist = sqrt(...
        periodic_diff(x1,x2,lx).^2 ...
      + periodic_diff(y1,y2,ly).^2 ...
      + periodic_diff(z1,z2,lz).^2);
end
function [dist] = calc_distance(c1, c2, lx, ly, lz)
  dist = sqrt(...
        periodic_diff(c1(1),c2(1),lx).^2 ...
      + periodic_diff(c1(2),c2(2),ly).^2 ...
      + periodic_diff(c1(3),c2(3),lz).^2);
end