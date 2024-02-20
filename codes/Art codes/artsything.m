close all; clc; clear;
%% Setting up random seed
rng(2023);

%% Setting up simulation parameters
ArtsyParams
dt = 0.05;
epsilon = quaternion(10^-8,0,0,0);
delta = 10^-8;

nprint = 10;
nsave = 20*5;
t0 = 5.0;%0.d0;
thold = 1000.d0;
savefigs = false;
savedata = true;

%% Generating grid points
[X, Y] = ndgrid(x0:dx:xf, y0:dy:yf);

%% Initializing order parameter n
n = zeros([Nx,Ny],'quaternion');
rand_pos = zeros(np, 2);
qlist = zeros(np,1,'quaternion');
nearest = zeros(Nx,Ny);

%%
curr_grain_id = 1;
exclude_distance = 25; %% 40 grid points

  rand_pos(:,1) = ceil(rand(np,1).*Nx);
  rand_pos(:,2) = ceil(rand(np,1).*Ny);
while curr_grain_id <= np
%   rand_num = rand();
%   temp = floor(Nx * Ny * rand_num);
%   rand_pos(curr_grain_id, 1) = floor(double(temp) / Nx) + 1;
%   rand_pos(curr_grain_id, 2) = mod(temp, Ny) + 1;
  D2=2;
  while D2>1
      a=2*rand()-1;
      b=2*rand()-1;
      c=2*rand()-1;
      d=2*rand()-1;
      D2= a*a+b*b+c*c+d*d;
  end
  qlist(curr_grain_id) = normalize(quaternion(a,b,c,d));
  curr_grain_id = curr_grain_id + 1;
end

for i = 1:Nx
    for j = 1:Ny
        dist_point = calc_distance(i,j,rand_pos(:,1), rand_pos(:,2), Nx, Ny);
        [~, nearest(i,j)] = min(dist_point);
    end
end
n = qlist(nearest);
figure(5) 
imagesc(nearest)
imagesc(norm(n-conj(n)))

%% Load
%load('n_ini.mat')
%load('n_ini_2.mat')
%load("n_pointy.mat")
load("art\archive\n_t=5.mat")
%load("art/n_t=235.mat")
%% init plots
figure(1)
til = tiledlayout(1,1);
til.Title.String = 't = 0';
nexttile(1)
h_dd = image(X(1, :), Y(:, 1), zeros(Nx,Ny,3));
hold on
h_mask = imshow(zeros(Nx,Ny,3), 'XData', X(1, :), 'YData', Y(:, 1));
hold off
axis image
set(gcf, 'position', [150, 100, 1400, 500])
%% Running simulations
t = t0;
steps = 0;
tic;
while t < thold 
  toc;
  fprintf('Initializing, t = %f\n', t)
  % update eta (n)
  dndt = zeros([Nx, Ny], 'quaternion');
  for istep = 1:nprint
    %% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    termn = 1-norm(n).^2;
    %normn = norm(n);

    lap_n = (circshift(n, [0, 1]) + circshift(n, [0, -1]) + ...
             circshift(n, [-1, 0]) + circshift(n, [1, 0]) - 4.d0 * n) ./ dx ./ dy;
    dndx = (circshift(n, [0, 1]) - circshift(n, [0, -1])) ./ (2*dx);
    dndy = (circshift(n, [1, 0]) - circshift(n, [-1, 0])) ./ (2*dy);
    %[dndx, dndy] = eno_gradArr_2D(n, termn, dx, dy);
    dndr2 = norm(dndx).^2 + norm(dndy).^2;
    
    terme  = (conj(n).*lap_n - conj(lap_n).*n); 
    
    %term1 = m0 * (-2.*(termn).*n);
    term1 = m0 * (-2.*(termn).*n);
    term1x = sqrt(1).*(conj(dndx).*term1 + conj(term1).*dndx).*dndx./norm(dndx).^2;
    term1x(isnan(term1x)) = quaternion(0,0,0,0);
    term1y = sqrt(1).*(conj(dndy).*term1 + conj(term1).*dndy).*dndy./norm(dndy).^2;
    term1y(isnan(term1y)) = quaternion(0,0,0,0);
    term1 = (term1x+term1y)./2;
    %term1 = m0 * (-2.*(termn).*n);
    term2 =  4.*(kappa/2).*(2.*n .*dndr2 - termn.*lap_n); %pointy one
    %term2 = (kappa/2).*(4.*n.*termn.*dndr2 - 2.* termn.^2 .*lap_n);
    %term2 = (kappa/2).*(4.*n.^2+2.*(1-n.^2)).*(-2).*((16.*n.^3+12.*termn.*n).*dndr2 + lap_n.*termn.*(4.*n.^2 + 2.*termn));
    %term2 = (kappa).*(2.*n.*dndr2);
    %term2 = 1.0*(kappa/2).*(-2.*lap_n).*termn;
    %term2 = (kappa/2).*(dndr2.*2.*n-2.*lap_n.*norm(n).^2);
    %term2 = (kappa/2).*(-2.*lap_n);

    dndt = -L * (2*term1 + 1*term2);
    %% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %% plot
    til.Title.String = ("t = " + num2str(t));
        nsq = norm(n).^2;
        C = reshape(compact(n),[Nx,Ny,4]);
        C = C(:,:,2:4);
    nexttile(1)
        max_mask = max(nsq,[],'all');
        min_mask = min(nsq,[],'all');
        temp_mask = 1 - (nsq - min_mask) ./ (max_mask - min_mask);
        set(h_dd, 'CData', (C+1)./2);
        set(h_mask, 'AlphaData', temp_mask)
    drawnow
    if savedata && mod(steps, nsave)==0
        save("art/n_t=" + num2str(t)+".mat", 'n');
    end
    if savefigs && mod(steps, nprint)==0
        print("../images/n_t=" + num2str(t), '-dpng', '-r500')
    end
    %% Update
    % update n and rho
    n = n + dndt * dt;
    
    % update time t
    t = t + dt;
    steps = steps+1;
  end 
end

%% Write n to files
save('../data/art_evolved.mat',"n");
fileID = fopen('../data/initial_n.dat','w');
fwrite(fileID, n, 'quaternion');
fclose(fileID);
%% defining functions
% function [dist] = calc_distance(x1, x2, y1, y2)
%   dist = sqrt( (x1 - x2)^2 + (y1 - y2)^2 );
% end

function [diff] = periodic_diff(x1, x2, lx)
  diff = (x1-x2);
  b = (abs(diff)>lx./2);
  diff = diff - sign(diff).*b.*(lx);
end
function [dist] = calc_distance(x1, y1, x2, y2, lx, ly)
  dist = sqrt(periodic_diff(x1,x2,lx).^2 + periodic_diff(y1,y2,ly).^2);
end

function C = qcolor(q)
    C = reshape(compact(q),[size(q),4]);
    C = C(:,:,2:4);
end