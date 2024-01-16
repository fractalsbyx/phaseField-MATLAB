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
[X, Y, Z] = meshgrid(x0:dx:xf, y0:dy:yf, z0:dz:zf);

%% Non-periodic Boundary
% bound_x=ones(Nx,Ny);
% bound_y=ones(Nx,Ny);
% bound_x([1 Nx],:)=0;
% bound_y(:,[1 Ny])=0;
% bound=bound_x.*bound_y;

%% Read in order parameters (grain information)
fileID = fopen('../data/initial_n.dat');
n = fread(fileID, Nx * Ny * Nz * np, 'double');
n = reshape(n, [Nx, Ny, Nz, np]);
fclose(fileID);

%% Mask Order Parameter
mask2 = ones(Nx, Ny, Nz);

%% Initialize Figures
if plotFigs
    figure(1)
    isodata = isosurface(X,Y,Z,n(:,:,:,1),0.5);
    iso = patch(isodata);
    %isonormals(X,Y,Z,n(:,:,:,1),iso)
    view(3);
    set(iso,'FaceColor',[0.5 1 0.5]);  
    set(iso,'EdgeColor','none');
    camlight;
    lighting gouraud;
    xlim([0,lx])
    ylim([0,ly])
    zlim([0,lz])
    drawnow
%     figure(2)
%     dummy = ones(Nx, Ny);
%     max_mask = 1;
%     min_mask = 0;
%     temp_mask = 1 - (dummy - min_mask) / (max_mask - min_mask);
%     
%     h_dd = image(X(1, :), Y(:, 1), rho, 'CDataMapping', 'scaled');
%     axis equal
%     black = zeros(Nx, Ny, 3);
%     hold on
%     h_mask = imshow(black, 'XData', X(1, :), 'YData', Y(:, 1));
%     hold off
%     set(h_mask, 'AlphaData', temp_mask)
%     xlim([x0 xf])
%     ylim([y0 yf])
%     xlabel('x^*')
%     ylabel('y^*')
%     title("title")
%     c = colorbar;
%     c.Label.String = 'Normalized Dislocation Density';
%     c.Label.FontSize = 24;
%     colormap(map)
%     set(gca, 'fontsize', 24)
%     set(gcf, 'position', [200 200 600 600])
%     clim([0*min_rho, max_rho])
%     axis on 
%     axis tight
%     drawnow
end

%% Starting simulation
sim_log = fopen("../data/sim_log.txt", 'w');
t = t0;
cycle = 0;
begin = tic;
sum_nsq = sum(n.^2, 4);
qx = quaternion(0,1,0,0);
qy = quaternion(0,0,1,0);
qz = quaternion(0,0,0,1);

while t <= tf
  fprintf('Simulating, t = %f\n', t);
  fprintf(sim_log, 'Simulating, t = %f\n', t);
  

  if saveFiles == true
    fileID = fopen("../data/n_t=" + num2str(t) + ".dat", 'w');
    fwrite(fileID, n, 'double');
    fclose(fileID);
%     fileID = fopen("../data/rho_t=" + num2str(t) + ".dat",'w');
%     fwrite(fileID, rho, 'double');
%     fclose(fileID);
  end
  
  if plotFigs
      [faces, verts] = isosurface(X,Y,Z,n(:,:,:,1),0.5);
      set(iso,'Faces',faces,'Vertices',verts);
%     % figure(2)
%     max_mask = max(sum_nsq(:));
%     min_mask = min(sum_nsq(:));
%     temp_mask = 1 - (sum_nsq - min_mask) / (max_mask - min_mask);
%     set(h_mask, 'AlphaData', temp_mask)
%     set(h_dd, 'CData', rho);
%     title(['t^* = ', num2str(t, '%.1f')])
  end
  
  if saveFigs && plotFigs
    print("../images/t=" + num2str(t), '-dpng', '-r500')
  end
  
  dndt = zeros(Nx, Ny, Nz, np);
  dndx = zeros(Nx, Ny, Nz, np);
  dndy = zeros(Nx, Ny, Nz, np);
  dndz = zeros(Nx, Ny, Nz, np);
  cycle = cycle + 1;
  tic;
  
  for istep = 1:nprint
    sum_nsq = sum(n.^2, 4);
%     [drhodx, drhody] = eno_gradArr_2D(rho, rho, dx, dy);
%     pairwise_sum = zeros(Nx, Ny);

    for i = 1:np
      [dndx(:,:,:,i), dndy(:,:,:,i), dndz(:,:,:,i)] = eno_gradArr_3D(n(:,:,:,i), n(:,:,:,i), dx, dy, dz);
    end
    qnormal = quaternion(zeros(Nx, Ny, Nz, np),dndx, dndy, dndz);

    for i = 1:np
      other_sum_nsq = sum_nsq - n(:,:,:,i).^2; 
      overlap_penalty = 2 * alpha * n(:,:,:,i) .* other_sum_nsq;
      % bulk term
      term1 = m0 * (-n(:,:,:,i) + n(:,:,:, i).^3 + overlap_penalty);
      
      % interfacial term
      lap_n = lap(n(:,:,:,i), dx,dy,dz);
      term2 = 0;
      term3 = -kappa *lap_n.*( ...
            norm(qnormal(:,:,:,i).*qx+conj(qx).*qnormal(:,:,:,i)) + ...
            norm(qnormal(:,:,:,i).*qy+conj(qy).*qnormal(:,:,:,i)) + ...
            norm(qnormal(:,:,:,i).*qz+conj(qz).*qnormal(:,:,:,i)) );
  
%       % stored energy term
%       grad_prod = (drhodx.*dndx(:,:,:,i) + drhody.*dndy(:,:,:,i));
%       term3 = beta * l_gb_ratio^2 * fs_const * grad_prod;

      % update dndt
      dndt(:,:,:,i) = -L * (term1 + term2 + term3);
    end
    %update n
    n = n + dndt.*dt;

    %update rho
%     term4 = (1/2).*(sum(abs(dndt), 3));
%     mask2 = mask2 - sign(mask2) .* term4 .* dt;
%     rho = (rho_init-rho_0) .* mask2.^2 ./ (mask2.^2 + (1-mask2).^2) + rho_0;

    % update time t
    t = t + dt;
  end
  tet = toc(begin);
  fprintf('Elapsed Time: %f\n', toc);
  fprintf(sim_log, 'Elapsed Time: %f\n', toc);
  fprintf('Total Elapsed Time: %f\n\n', tet);
  fprintf(sim_log, 'Total Elapsed Time: %f\n\n', tet);
end
fclose(sim_log);
