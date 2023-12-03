close all; clc; clear;
%% Setting up random seed
rng(2022);

%% Setting up simulation parameters
Params
dt = 0.1;

nprint = 100;
t0 = 0.d0;
thold = 50.d0;

%% Generating grid points
[X, Y] = meshgrid(y0:dy:yf, x0:dx:xf);

%% Initializing order parameter n
n = zeros(Nx, Ny, np);
range = 3; % 3
rand_pos = zeros(np, 2);

curr_grain_id = 1;
exclude_distance = 25; %% 40 grid points
while curr_grain_id <= np
  rand_num = rand();
  temp = floor(Nx * Ny * rand_num);
  rand_pos(curr_grain_id, 1) = floor(double(temp) / Ny) + 1;
  rand_pos(curr_grain_id, 2) = mod(temp, Nx) + 1;
  
  is_valid_grain = true;
  for past_grain_id = 1 : (curr_grain_id - 1)
    min_x = min(rand_pos(curr_grain_id, 1), rand_pos(past_grain_id, 1));
    max_x = max(rand_pos(curr_grain_id, 1), rand_pos(past_grain_id, 1));
    min_y = min(rand_pos(curr_grain_id, 2), rand_pos(past_grain_id, 2));
    max_y = max(rand_pos(curr_grain_id, 2), rand_pos(past_grain_id, 2));
    
    dist1 = calc_distance(min_x, max_x, min_y, max_y);
    dist2 = calc_distance(min_x, max_x - Nx, min_y, max_y);
    dist3 = calc_distance(min_x, max_x, min_y, max_y - Ny);
    dist4 = calc_distance(min_x, max_x - Nx, min_y, max_y - Ny);
    
    distance = min([dist1, dist2, dist3, dist4]);
    
    if distance < exclude_distance
      fprintf("curr pos = [%i, %i], past pos = [%i, %i], distance = %f\n", rand_pos(curr_grain_id, 1), rand_pos(curr_grain_id, 2), rand_pos(past_grain_id, 1), rand_pos(past_grain_id, 2), distance);
      is_valid_grain = false;
      break;
    end
  end
  
  if is_valid_grain == true
    for x_id = (rand_pos(curr_grain_id, 1) - range) : (rand_pos(curr_grain_id, 1) + range)
      for y_id = (rand_pos(curr_grain_id, 2) - range) : (rand_pos(curr_grain_id, 2) + range)
        if calc_distance(x_id, rand_pos(curr_grain_id, 1), y_id, rand_pos(curr_grain_id, 2)) <= range
          curr_x = 1+(mod(x_id, Nx));
          curr_y = 1+(mod(y_id, Ny));
%           if x_id < 1
%             curr_x = x_id + Nx;
%           elseif x_id > Nx
%             curr_x = x_id - Nx;
%           else
%             curr_x = x_id;
%           end
%           if y_id < 1
%             curr_y = y_id + Ny;
%           elseif y_id > Ny
%             curr_y = y_id - Ny;
%           else
%             curr_y = y_id;
%           end
          n(curr_x, curr_y, curr_grain_id) = 1.0;
        end
      end
    end
    
    curr_grain_id = curr_grain_id + 1;
  end
  
end

%% Running simulations
t = t0;
tic;
while t < thold 
  toc;
  fprintf('Initializing, t = %f\n', t)
  
  % update eta (n)
  dndt = zeros(Nx, Ny, np);
  for istep = 1:nprint
    all_sum_nsq = sum(n.^2, 3);

    for i = 1:np
      sum_nsq = all_sum_nsq - n(:, :, i).^2; 
      term_sum_nsq = 2 * alpha * n(:, :, i) .* sum_nsq;
      % bulk term
      term1 = m0 * (-n(:, :, i) + n(:, :, i).^3 + term_sum_nsq);
      
      % interfacial term
      lap_n = (circshift(n(:, :, i), [0, 1]) + circshift(n(:, :, i), [0, -1]) + ...
               circshift(n(:, :, i), [-1, 0]) + circshift(n(:, :, i), [1, 0]) - 4.d0 * n(:, :, i)) / dx / dy;
      term2 = -kappa * lap_n;

      % update dndt
      dndt(:, :, i) = -L * (term1 + term2);
    end
    
    % update n and rho
    n = n + dndt * dt;
    
    % update time t
    t = t + dt;
  end 
end

% %% Initializaing dislocation density rho
% all_sum_nsq = sum(n.^2, 3);
% rho = zeros(Nx, Ny);
% for i = 1:length(rho_const)
%   rho = rho + n(:, :, i).^2 .* rho_const(i);
% end
% rho = rho ./ all_sum_nsq;

%% Write n to files
fileID = fopen('../data/initial_n.dat','w');
fwrite(fileID, n, 'double');
fclose(fileID);
fprintf("Initial condition written to ../data/initial_n.dat\n")

% %% Write rho to files
% fileID = fopen('../data/initial_rho.dat','w');
% fwrite(fileID, rho, 'double');
% fclose(fileID);
% 
% fileID = fopen('../data/rho_const.dat','w');
% fwrite(fileID, rho_const, 'double');
% fclose(fileID);

%% defining functions
function [dist] = calc_distance(x1, x2, y1, y2)
  dist = sqrt( (x1 - x2)^2 + (y1 - y2)^2 );
end