close all; clc; clear;
%% Setting up random seed
rng(42);

map = jet;
map = map(70:220, :);

%% Setting up simulation parameters
Params

t0 = 0;
step = 1;
tf = 10;


%% Generating grid points
[X, Y] = meshgrid(y0:dy:yf, x0:dx:xf);

%% Plotting figures
for t = t0:step:tf
  
  fileID = fopen("../data/n_t=" + num2str(t) + ".dat");
  n = fread(fileID, Nx * Ny * np, 'double');
  n = reshape(n, [Nx, Ny]);
  fclose(fileID);

  h_dd = image(X(1, :), Y(:, 1), n, 'CDataMapping', 'scaled');
  axis square
  %close all
  xlim([x0 xf])
  ylim([y0 yf])
  xlabel('x^*')
  ylabel('y^*')
  title(['t^* = ', num2str(t)])
  c = colorbar;
  c.Label.String = 'Orientation';
  c.Label.FontSize = 24; % 24
  colormap(map)
  set(gca, 'fontsize', 24) % 24
  set(gcf, 'position', [200 200 700 500])
  % set(gcf, 'position', [200 200 1000 650])
  % clim([0*min_rho, max_rho])
  axis on
  axis tight
  drawnow
  
  print("../images/t=" + num2str(t), '-dpng', '-r300')
end









