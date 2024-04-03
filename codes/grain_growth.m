close all; clc; clear;
%% Setting up random seed
rng(2021);

%% Version of Driving Force
mode = 'gradprod';

map = jet;
map = map(70:220, :);

%% Setting up simulation parameters
Params % See Params.m
dt = 0.02; % 0.1

nprint = 10; % 50; % 1/dt;
t0 = 0.d0;
tf = 10.0 +dt; % +dt so last frame will print

saveFiles = true;
saveFigs = true;
plotFigs = true;

%% Generating grid points
X = x0:dx:xf;

%% Non-periodic Boundary
bound_x=ones(Nx,1);
bound_x([1 Nx],:)=0;

%% Read in order parameters (grain information)
fileID = fopen('../data/initial_n.dat');
n = fread(fileID, Nx * np, 'double');
n = reshape(n, [Nx, np]);
fclose(fileID);

%% Mask Order Parameter
mask2 = ones(Nx, 1);

%% Read in dislocation density field (rho)
rho_grain_id
fileID = fopen('../data/rho_const.dat'); % Gentry Model
rho_const = fread(fileID, np, 'double');
fclose(fileID);

fileID = fopen('../data/initial_rho.dat'); % Gentry Model
rho = fread(fileID, Nx, 'double');
fclose(fileID);

%% Normalize Rho
rho_ref = 10^15;%mean(rho, "all");
fs_const = 0.5 * G * b^2 * rho_ref / m0_const;
fs_const = 0.1;

min_rho = min(rho_const); % for figures
max_rho = max(rho_const);

%% Initialize Figures
if plotFigs
    figure(2)
    hold on
    plt1 = plot(X, n(:,1), '-o');
    plt2 = plot(X, n(:,2), '-o');
    pltrho = plot(X,rho);
    axis on
    axis tight
    drawnow
end

%% Starting simulation
t = t0;
cycle = 1;
begin = tic;
all_sum_nsq = sum(n.^2, 3);
vrecord = zeros(floor((tf-t0)/dt) + 1,1);

while t <= tf
    fprintf('Simulating, t = %f\n', t);

    if saveFiles == true
        fileID = fopen("../data/n_t=" + num2str(t) + ".dat", 'w');
        fwrite(fileID, n, 'double');
        fclose(fileID);

        fileID = fopen("../data/rho_t=" + num2str(t) + ".dat",'w');
        fwrite(fileID, rho, 'double');
        fclose(fileID);
    end

    if plotFigs
        set(plt1, 'YData', n(:,1));
        set(plt2, 'YData', n(:,2));
        set(pltrho, 'YData', rho);
        title(['t^* = ', num2str(t, '%.1f')])
    end

    if plotFigs && saveFigs
        print("../images/t=" + num2str(t), '-dpng', '-r500')
    end

    dndt = zeros(Nx, np);
    dndx = zeros(Nx, np);

    tic;

    for istep = 1:nprint
        nsq = n.^2;
        all_sum_nsq = sum(n.^2, 2);
        all_sq_sum_nsq = all_sum_nsq.^2;
        drhodx = gradArr_1D(rho, rho, dx);

        for i = 1:np
            dndx(:, i) = gradArr_1D(n(:, i), rho, dx).*bound_x;
        end
        termx = zeros(Nx, np);
        for i = 1:np
            sum_nsq = all_sum_nsq - n(:, i).^2;
            term_sum_nsq = 2 * alpha * n(:, i) .* sum_nsq;
            % bulk term
            term1 = m0 * (-n(:, i) + n(:, i).^3 + term_sum_nsq);

            % interfacial term
            lap_n = lap_1D(n(:,i), dx).*bound_x;
            term2 = -kappa * lap_n;

            % stored energy term
            if     strcmp(mode, 'gradprod')
                grad_prod = drhodx.*dndx(:, i);
                term3 = beta * l_gb_ratio^2 * fs_const * grad_prod;
            elseif strcmp(mode, 'gentry')
                term3 = fs_const * 2*n(:,i)./all_sum_nsq.^2 .* (rho_const(i)-rho);
            elseif strcmp(mode, 'cubic')
                term3 = fs_const * (6*n(:,i) - 6*n(:,i).^2) .* (rho_const(i)-rho);
            end


        % update dndt
        dndt(:, i) = -L * (term1 + term2 + term3);

        if i==1
            V1 = dx*sum(dndt(:,1));
            vrecord(cycle) = V1;
        end
        dndt(:,i) = dndt(:,i) + V1*dndx(:, i);


        end
        %update n
        n = n + dndt.*dt;

        %update rho
        term5 = (1/2).*(sum(abs(dndt), 3));
        mask2 = mask2 - sign(mask2) .* term5 .* dt;
        rho = zeros(Nx, 1);
        for i = 1:np
            rho = rho + rho_const(i)*n(:,i).*n(:,i);
        end
        rho = rho./all_sum_nsq;

        % update time t
        t = t + dt;
        cycle = cycle + 1;
    end
    tet = toc(begin);
    fprintf('Elapsed Time: %f\n', toc);
    fprintf('Total Elapsed Time: %f\n\n', tet);
end
vrecord = vrecord(1:500);
fileID = fopen("../data/IntegratedForce_" + mode + ".dat", 'w');
fwrite(fileID, vrecord, 'double');
fclose(fileID);