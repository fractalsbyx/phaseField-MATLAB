close all; clc; clear;
%% Setting up random seed
rng(2023);

%% Setting up simulation parameters
Params
dt = 0.025;
epsilon = quaternion(10^-8,0,0,0);
delta = 10^-8;

nprint = 100;
t0 = 0.d0;
thold = 1000.d0;

%% Generating grid points
X = x0:dx:xf;

%% Initializing order parameter n
n = zeros([Nx,1],'quaternion');
range = 9; % 3
rand_pos = zeros(np, 2);

%%
r1 = normalize(quaternion(rand-0.5,rand-0.5,rand-0.5,rand-0.5));
r2 = normalize(quaternion(rand-0.5,rand-0.5,rand-0.5,rand-0.5));
n(X>Nx*dx/2) = r1;%quaternion(1,0,0,0);
n(X<Nx*dx/3) = r2;quaternion(-1,0,0,0);
%n = (1/2).*(n-conj(n));

%% Load
load('n_ini.mat')
n = n(100,:);
%load('n_ini_2.mat')
%% init plots
figure(1)
til = tiledlayout(1,2);
til.Title.String = 't = 0';
nexttile(1)
h_dd = plot(X, zeros(Nx,1));
ylim([-1.25,1.25])
axis square
nexttile(2)
[sX,sY,sZ] = sphere;
%web = scatter3(zeros(Nx*Ny,1),zeros(Nx*Ny,1),zeros(Nx*Ny,1),1,zeros(Nx*Ny,3));
web2 = plot3(zeros(Nx,1),zeros(Nx,1),zeros(Nx,1));
hold on
quiv = quiver3(zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny),zeros(Nx,Ny));
mesh(sX,sY,sZ,'FaceAlpha',0.1,'EdgeAlpha',0.1,'FaceColor','r');
hold off
colorbar
xlim([-1,1]);
ylim([-1,1]);
zlim([-1,1]);
axis square
set(gcf, 'position', [150, 100, 1000, 500])
%% Running simulations
t = t0;
tic;
while t < thold 
  toc;
  fprintf('Initializing, t = %f\n', t)
  % update eta (n)
  dndt = zeros([Nx, 1], 'quaternion');
  for istep = 1:nprint
    %% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    termn = 1-norm(n).^2;
    n_normal = normalize(n);
    n_normal(isnan(n_normal)) = quaternion(0,0,0,0);

    lap_n = (circshift(n, 1) + circshift(n, -1) - 2.d0 * n) ./ dx ./dx;
    dndx  = (circshift(n, 1) - circshift(n, -1)) ./ (2*dx);
    [dndx, ~]  = eno_gradArr_2D(n,-termn,dx,1);
    dndx_normal = normalize(n);
    dndx_normal(isnan(dndx_normal)) = quaternion(0,0,0,0);
    dndr2 = norm(dndx).^2;
    normality = kappa.*reshape(compact((dndx./n).^2),Nx,4);
    normality = normality(:,1);
    [ax,bx,cx,fx] = parts(n./dndx.^2);

    termax = conj(n).*dndx - conj(dndx).*n;
    termbx = dndx - conj(dndx);
    termcx = 2*(termbx.*termax + termax.*termbx);
    termd  = conj(n)-n;
    terme  = (conj(n).*lap_n - conj(lap_n).*n); 
    termf  = termd.*terme + terme.*termd;
    termg  = (-1/4).*(termcx-termf);
    termh  = 4.*termn .*termg + (1/4).*(norm(termax).^2).*(16).*n;

    term3 = (1*kappa/2).*(termg);

    %term1 = m0 * (-2.*(termn).*n);
    n_parallel = normalize(conj(dndx)).*n;  %adj 3
    c_line = center(n_parallel);                 
    n_shifted = -c_line;    %adj 1 & 2
    term1_shifted = m0 * (-2.*(1-norm(n_shifted).^2).*n_shifted);
    term1_adjustment_1 = mobius_v(-n_shifted, n_shifted, term1_shifted);
    term1_adjustment_2 = mobius_v(n_parallel, quaternion(0,0,0,0), term1_adjustment_1);
    term1_adjustment_3 = normalize(dndx).* term1_adjustment_2;
    term1 = term1_adjustment_3;

    lap_n_shifted = (circshift(n_shifted, 1) + circshift(n_shifted, -1) - 2.d0 * n_shifted) ./ dx ./dx;
    term2_shifted = (kappa/2).*(-2.*lap_n);
    term2_adjustment_1 = mobius_v(-n_shifted, n_shifted, term1_shifted);
    term2_adjustment_2 = mobius_v(n_parallel, quaternion(0,0,0,0), term2_adjustment_1);
    term2_adjustment_3 = normalize(dndx).* term2_adjustment_2;
    term2 = term2_adjustment_3;

%     term1a = m0 * (-2.*termn.*n);
%     term1 = (conj(dndx).*term1a + conj(term1a).*dndx).*dndx./norm(dndx).^2;
     term1(isnan(term1)) = quaternion(0,0,0,0);
     term2(isnan(term2)) = quaternion(0,0,0,0);
     term1 = m0 * Dot(n_normal, dndx_normal) .*(-2.*(termn.^2).*n);
    %term1 = 1*m0 * (-2.*(termn.^2).*n);
    %term2 = (kappa/2).*(2.*n .*dndr2 - termn.*lap_n);
    %term2 = (kappa/2).*(4.*n.*termn.*dndr2 - 2.* termn.^2 .*lap_n);
    %term2 = (kappa/2).*(4.*n.^2+2.*(1-n.^2)).*(-2).*((16.*n.^3+12.*termn.*n).*dndr2 + lap_n.*termn.*(4.*n.^2 + 2.*termn));
    %term2 = (kappa).*(2.*n.*dndr2);
    term2 = Dot(n_normal, dndx_normal) .*(kappa/2).*(-2.*lap_n);
    %term2 = (kappa/2).*(dndr2.*2.*n-2.*lap_n.*norm(n).^2);

    dndt = -L * (1*term1 + 1*term2 + 0*term3+1e-10);
    %% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %% plot
    til.Title.String = ("t = " + num2str(t));
        nqt = norm(n).^4;
        myv = kappa.*(norm(dndx.*conj(n)-conj(n).*dndx).^2);
        V = reshape(compact(terme),[Nx,4]);
        V = V(:,2:4);
        C = reshape(compact(n),[Nx,4]);
        C0 = C(:,1);
        C = C(:,2:4);
    nexttile(1)
        %set(h_dd, 'CData', (C+1)./2);
        set(h_dd, 'YData', C0);
    nexttile(2)
        C2 = (reshape(C,Nx,3));
        set(web2, 'XData',C(:,1),'YData',C(:,2),'ZData',C(:,3));
        %set(quiv, 'XData',C(:,1),'YData',C(:,2),'ZData',C(:,3),'UData',V(:,1),'VData',V(:,2),'WData',V(:,3))
    drawnow
    %% Update
    % update n and rho
    n = n + dndt * dt;
    
    % update time t
    t = t + dt;
  end 
end

%% Write n to files
save('../data/n_evolved.mat',"n");
fileID = fopen('../data/initial_n.dat','w');
fwrite(fileID, n, 'quaternion');
fclose(fileID);
%% defining functions
function [dist] = calc_distance(x1, x2, y1, y2)
  dist = sqrt( (x1 - x2)^2 + (y1 - y2)^2 );
end

function y = mobius(shift, point)
    y = (point + shift)./(1+shift.*conj(point));
end

function y = mobius_v(shift, point, vector)
    y = ((1-norm(shift).^2)./(1+shift.*point).^2) .* vector;
end

function y = center(point)
    S = norm(point);
    r = point + conj(point);
    a = r.*(S.*S - 1);
    b = (S.^4 - 1);
    y = (-b - sqrt(norm(b.*b - 4.*a.*a)))./(2.*a);
end

function y = Dot(a, b)
    y = conj(b).*a + conj(a).*b;
end