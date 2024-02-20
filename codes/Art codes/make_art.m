close all; clc; clear;
%% Setting up random seed
rng(42);

map = jet;
map = map(70:220, :);

%% Setting up simulation parameters
ArtsyParams

t0 = 5;
step = 5;
tf = 60;


%% Generating grid points
[X, Y] = meshgrid(y0:dy:yf, x0:dx:xf);

%% Plotting figures
for t = t0:step:tf
  
%   fileID = fopen("../data/n_t=" + num2str(t) + ".dat");
%   n = fread(fileID, Nx * Ny * np, 'double');
%   n = reshape(n, [Nx, Ny, np]);
%   fclose(fileID);
  load("art\n_t="+num2str(t));
  imwrite(qcolor(n), "art\n_t="+num2str(t)+".png");
%   all_sum_nsq = sum(n.^2, 3);
  
%   max_mask = max(all_sum_nsq(:));
%   min_mask = min(all_sum_nsq(:));
%   temp_mask = 1 - (all_sum_nsq - min_mask) / (max_mask - min_mask);
%   
%   %pcolor(temp_mask)
% 
%   h_dd = image(X(1, :), Y(:, 1), rho, 'CDataMapping', 'scaled');
%   axis square
%   %close all
%   black = cat(3, zeros(size(all_sum_nsq)), zeros(size(all_sum_nsq)), zeros(size(all_sum_nsq)));
%   hold on
%   h_mask = imshow(black, 'XData', X(1, :), 'YData', Y(:, 1));
%   hold off
%   set(h_mask, 'AlphaData', temp_mask)
%   xlim([x0 xf])
%   ylim([y0 yf])
%   xlabel('x^*')
%   ylabel('y^*')
%   title(['t^* = ', num2str(t)])
%   c = colorbar;
%   c.Label.String = 'Normalized Dislocation Density';
%   c.Label.FontSize = 24; % 24
%   colormap(map)
%   set(gca, 'fontsize', 24) % 24
%   set(gcf, 'position', [200 200 700 500])
% %   set(gcf, 'position', [200 200 1000 650])
%   clim([0*min_rho, max_rho])
%   axis on
%   axis tight
%   drawnow
%   
%   print("../images/t=" + num2str(t), '-dpng', '-r300')
end
function C = qcolor(q)
    C = reshape(compact(q),[size(q),4]);
    C = (C(:,:,2:4)+1)./2;
end