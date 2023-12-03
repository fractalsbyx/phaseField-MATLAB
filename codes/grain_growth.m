close all; clc; clear;
%% Setting up random seed
rng(2021);

map = jet;
map = map(70:220, :);

%% Setting up simulation parameters
Params % See Params.m
dt = 0.1; % 0.1

nprint = 10; % 50; % 1/dt;
t0 = 0.d0;
tf = 10.0 +1; % +1 so last frame will print

saveFiles = true;
saveFigs = true;
plotFigs = true;

%% Generating grid points
[X, Y] = meshgrid(y0:dy:yf, x0:dx:xf);

%% Non-periodic Boundary
% bound_x=ones(Nx,Ny);
% bound_y=ones(Nx,Ny);
% bound_x([1 Nx],:)=0;
% bound_y(:,[1 Ny])=0;
% bound=bound_x.*bound_y;

%% Read in order parameters (grain information)
fileID = fopen('../data/initial_n.dat');
n = fread(fileID, Nx * Ny * np, 'double');
n = reshape(n, [Nx, Ny, np]);
fclose(fileID);

%% Mask Order Parameter
mask2 = ones(Nx, Ny);

%% Read in dislocation density field (rho)

fileID = fopen('../data/initial_rho.dat');
rho = fread(fileID, Nx*Ny, 'double');
rho = reshape(rho, [Nx, Ny]);
fclose(fileID);

%% Proof-of-Concept Things
% rho_linear
% sin_square
% sin_sum

% fileID = fopen('../data/rho_linear.dat'); % linear gradient
% rho = fread(fileID, Nx*Ny, 'double');
% rho = reshape(rho, [Nx, Ny]);
% fclose(fileID);

% fileID = fopen('../data/rho_interp.dat'); % gradient from experiment
% rho = fread(fileID, Nx*Ny, 'double');
% rho = reshape(rho, [Nx, Ny]);
% fclose(fileID);

% fileID = fopen('../data/data_smoothnoise.dat'); % just smoothed noise
% smooth = fread(fileID, Nx*Ny, 'double');
% smooth = reshape(smooth, [Nx, Ny]);
% fclose(fileID);

% fileID = fopen("../data/sin_sum.dat"); % smooth localized shape
% sin_sqr = fread(fileID, Nx*Ny, 'double');
% sin_sqr = reshape(sin_sqr, [Nx, Ny]);
% fclose(fileID);

% rho = sin_sqr;
% rho_mean = mean(rho(:));
% rho = rho_mean + rho_mean*(0.2*smooth+1*sin_sqr);

%% Normalize Rho
rho_mean = mean(rho, "all");
fs_const = 0.5 * G * b^2 * rho_mean / m0_const;

fid_energy = fopen("../data/sim_log.txt", 'w');
fprintf("rho_mean = %i\n\n", rho_mean);
fprintf(fid_energy, "rho_mean = %i\n\n", rho_mean);

rho = rho / rho_mean; % non-dim
rho_0 = rho_0 / rho_mean;
rho_init = rho;

min_rho = min(rho,[],'all'); % for figures
max_rho = max(rho,[],'all');

% fileID = fopen('../data/rho_const.dat'); % Gentry Model 
% rho_const = fread(fileID, np, 'double');
% fclose(fileID);
% rho_const = rho_const / rho_mean; % non-dim

%% Initialize Figures
if plotFigs
    figure(2)
    dummy = ones(Nx, Ny);
    max_mask = 1;
    min_mask = 0;
    temp_mask = 1 - (dummy - min_mask) / (max_mask - min_mask);
    
    h_dd = image(X(1, :), Y(:, 1), rho, 'CDataMapping', 'scaled');
    axis equal
    black = zeros(Nx, Ny, 3);
    hold on
    h_mask = imshow(black, 'XData', X(1, :), 'YData', Y(:, 1));
    hold off
    set(h_mask, 'AlphaData', temp_mask)
    xlim([x0 xf])
    ylim([y0 yf])
    xlabel('x^*')
    ylabel('y^*')
    title("title")
    c = colorbar;
    c.Label.String = 'Normalized Dislocation Density';
    c.Label.FontSize = 24;
    colormap(map)
    set(gca, 'fontsize', 24)
    set(gcf, 'position', [200 200 600 600])
    clim([0*min_rho, max_rho])
    axis on 
    axis tight
    drawnow
end

%% Starting simulation
t = t0;
cycle = 0;
begin = tic;
all_sum_nsq = sum(n.^2, 3);

while t <= tf
  fprintf('Simulating, t = %f\n', t);
  fprintf(fid_energy, 'Simulating, t = %f\n', t);
  

  if saveFiles == true
    fileID = fopen("../data/n_t=" + num2str(t) + ".dat", 'w');
    fwrite(fileID, n, 'double');
    fclose(fileID);

    fileID = fopen("../data/rho_t=" + num2str(t) + ".dat",'w');
    fwrite(fileID, rho, 'double');
    fclose(fileID);
  end
  
  if plotFigs
    % figure(2)
    max_mask = max(all_sum_nsq(:));
    min_mask = min(all_sum_nsq(:));
    temp_mask = 1 - (all_sum_nsq - min_mask) / (max_mask - min_mask);
    set(h_mask, 'AlphaData', temp_mask)
    set(h_dd, 'CData', rho);
    title(['t^* = ', num2str(t, '%.1f')])
  end
  
  if saveFigs && plotFigs
    print("../images/t=" + num2str(t), '-dpng', '-r500')
  end
  
  dndt = zeros(Nx, Ny, np);
  dndx = zeros(Nx, Ny, np);
  dndy = zeros(Nx, Ny, np);
  cycle = cycle + 1;
  tic;
  
  for istep = 1:nprint
    nsq = n.^2;
    all_sum_nsq = sum(n.^2, 3);
    all_sq_sum_nsq = all_sum_nsq.^2;
    [drhodx, drhody] = eno_gradArr_2D(rho, rho, dx, dy);
    pairwise_sum = zeros(Nx, Ny);

    for i = 1:np
      [dndx(:, :, i), dndy(:, :, i)] = eno_gradArr_2D(n(:, :, i), rho, dx, dy);
    end

    termx = zeros(Nx, Ny, np);
    termy = zeros(Nx, Ny, np);
    for i = 1:np
      sum_nsq = all_sum_nsq - n(:, :, i).^2; 
      term_sum_nsq = 2 * alpha * n(:, :, i) .* sum_nsq;
      % bulk term
      term1 = m0 * (-n(:, :, i) + n(:, :, i).^3 + term_sum_nsq);
      
      % interfacial term
      lap_n = (circshift(n(:, :, i), [0, 1]) + circshift(n(:, :, i), [0, -1]) + ...
               circshift(n(:, :, i), [-1, 0]) + circshift(n(:, :, i), [1, 0]) - 4.d0 * n(:, :, i)) / dx / dy;
      term2 = -kappa * lap_n;
  
      % stored energy term
      grad_prod = (drhodx.*dndx(:, :, i) + drhody.*dndy(:, :, i));
      term3 = beta * l_gb_ratio^2 * fs_const * grad_prod;

      % update dndt
      dndt(:, :, i) = -L * (term1 + term2 + term3);
    end
    %update n
    n = n + dndt.*dt;

    %update rho
    term4 = (1/2).*(sum(abs(dndt), 3));
    mask2 = mask2 - sign(mask2) .* term4 .* dt;
    rho = (rho_init-rho_0) .* mask2.^2 ./ (mask2.^2 + (1-mask2).^2) + rho_0;

    % update time t
    t = t + dt;
  end
  tet = toc(begin);
  fprintf('Elapsed Time: %f\n', toc);
  fprintf(fid_energy, 'Elapsed Time: %f\n', toc);
  fprintf('Total Elapsed Time: %f\n\n', tet);
  fprintf(fid_energy, 'Total Elapsed Time: %f\n\n', tet);
end
fclose(fid_energy);


%% old code
%     term5 = zeros(Nx, 1);
%     for i=1:np
%         term5 = term5 + n(:,:, i) .* (1-n(:,:, i));
%     end

%     drho_maskdt = -2.*term4.*pairwise_sum.*sign(rho_mask)./all_sq_sum_nsq;
%     drho_maskdt =    -term4.*term5       .*sign(rho_mask)./all_sq_sum_nsq;
%     rho_mask = rho_mask + drho_maskdt.*dt;
%     rho = (rho_init-rho_0) .* rho_mask+rho_0;

