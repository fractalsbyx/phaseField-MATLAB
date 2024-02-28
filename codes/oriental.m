close all; clc; clear;
%% Setting up random seed
rng(2022);

%% Setting up simulation parameters
Params
dt = 0.05;
epsilon = quaternion(10^-8,0,0,0);
delta = 10^-8;

nprint = 10;
t0 = 0.d0;
thold = 1000.d0;
savefigs = false;
savedata = false;

%% Generating grid points
[X, Y] = meshgrid(y0:dy:yf, x0:dx:xf);

%% Load
% load('n_ini.mat')
% load('n_ini_2.mat')
% load("n_pointy.mat")
load("../data/initial_n.mat")
%% init plots
figure(1)
til = tiledlayout(1,3);
til.Title.String = 't = 0';
nexttile(1)
h_dd = image(X(1, :), Y(:, 1), zeros(Nx,Ny,3));
hold on
h_mask = imshow(zeros(Nx,Ny,3), 'XData', X(1, :), 'YData', Y(:, 1));
hold off
axis square
nexttile(2)
mag = imagesc(zeros(Nx,Ny));
colorbar
colormap('turbo');
axis square
nexttile(3)
%web = scatter3(zeros(Nx*Ny,1),zeros(Nx*Ny,1),zeros(Nx*Ny,1),1,zeros(Nx*Ny,3));
web2 = mesh(zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny,3),'FaceAlpha',0.0);
hold on
quiv = quiver3(zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny));
[sX,sY,sZ] = sphere;
mesh(sX,sY,sZ,'FaceAlpha',0.1,'EdgeAlpha',0.1,'FaceColor','c');
hold off
%set(p,'FaceAlpha',0.4)
colorbar
xlim([-1,1]);
ylim([-1,1]);
zlim([-1,1]);
axis square
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
    normn = norm(n);

    lap_n = (circshift(n, [0, 1]) + circshift(n, [0, -1]) + ...
             circshift(n, [-1, 0]) + circshift(n, [1, 0]) - 4.d0 * n) ./ dx ./ dy;
    dndx = (circshift(n, [0, 1]) - circshift(n, [0, -1])) ./ (2*dx);
    dndy = (circshift(n, [1, 0]) - circshift(n, [-1, 0])) ./ (2*dy);
    %[dndx, dndy] = eno_gradArr_2D(n, -normn, dx, dy);
    dndr2 = norm(dndx).^2 + norm(dndy).^2;

    termax = conj(n).*dndx - conj(dndx).*n;
    termay = conj(n).*dndy - conj(dndy).*n;
    termbx = dndx - conj(dndx);
    termby = dndy - conj(dndy);
    termcx = 2*(termbx.*termax + termax.*termbx);
    termcy = 2*(termby.*termay + termay.*termby);
    termd  = conj(n)-n;
    terme  = (conj(n).*lap_n - conj(lap_n).*n); 
    termf  = termd.*terme + terme.*termd;
    termg  = (-1/4).*(termcx+termcy-termf);
    termh  = 4.*termn.*termg + (1/4).*(norm(termax).^2 + norm(termay).^2).*(16).*n;
    %normality = kappa.*reshape(compact((dndx./n).^2 +(dndy./n).^2),Nx,Ny,4);
    %normality = normality(:,:,1);
    normality = kappa.*norm((termax).^2 + (termay).^2);

    term3 = (1*kappa/2).*(termg);

    %term1 = m0 * (-2.*(termn).*n);
    term1a = m0 * (-2.*(termn).*n);
    term1x = sqrt(1).*(conj(dndx).*term1a + conj(term1a).*dndx).*dndx./norm(dndx).^2;
    term1x(isnan(term1x)) = quaternion(0,0,0,0);
    term1y = sqrt(1).*(conj(dndy).*term1a + conj(term1a).*dndy).*dndy./norm(dndy).^2;
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

    dndt = -L * (2*term1 + 1*term2 + 0*term3);
    %% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %% plot
    til.Title.String = ("t = " + num2str(t));
        nsq = norm(n).^2;
        myv = kappa.*(norm(dndx.*conj(n)-conj(n).*dndx).^2 + norm(dndy.*conj(n)-conj(n).*dndy).^2);
        V = reshape(compact(terme),[Nx,Ny,4]);
        V = V(:,:,2:4);
        C = reshape(compact(n),[Nx,Ny,4]);
        C0 = C(:,:,1);
        C = C(:,:,2:4);
    nexttile(1)
        max_mask = max(nsq,[],'all');
        min_mask = min(nsq,[],'all');
        temp_mask = 1 - (nsq - min_mask) ./ (max_mask - min_mask);
        set(h_dd, 'CData', (C+1)./2);
        set(h_mask, 'AlphaData', temp_mask)
    nexttile(2)
        set(mag, 'CData', nsq)%nsq
    nexttile(3);
        C2 = (reshape(C,Nx*Ny,3));
        %set(web, 'XData',C2(:,1),'YData',C2(:,2),'ZData',C2(:,3),'CData',dndr2(:));
        %set(web2, 'XData',C(:,:,1),'YData',C(:,:,2),'ZData',C(:,:,3),'CData',C0);
        set(web2, 'XData',C(:,:,1),'YData',C(:,:,2),'ZData',C(:,:,3),'CData',myv);
        set(quiv, 'XData',C(:,:,1),'YData',C(:,:,2),'ZData',C(:,:,3),'UData',V(:,:,1),'VData',V(:,:,2),'WData',V(:,:,3))
    drawnow
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
save('../data/n_evolved.mat',"n");