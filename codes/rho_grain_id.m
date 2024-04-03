%% Setting up random seed
rng(2021);

%% Params
Params

%% Read in order parameters (grain information)
fileID = fopen('../data/initial_n.dat');
n = fread(fileID, Nx * np, 'double');
n = reshape(n, [Nx, np]);
fclose(fileID);

%% Rho Constants
rho_const = normrnd(rho_mean, rho_std, np, 1);
rho_const(1) = 0;
rho_const(2) = 1;
%% Initializaing dislocation density rho
all_sum_nsq = sum(n.^2, 2);
rho = zeros(Nx, 1);
for i = 1:length(rho_const)
  rho = rho + n(:, i).^2 .* rho_const(i);
end
rho = rho ./ all_sum_nsq;

%% Write rho to files
fileID = fopen('../data/initial_rho.dat','w');
fwrite(fileID, rho, 'double');
fclose(fileID);

fileID = fopen('../data/rho_const.dat','w');
fwrite(fileID, rho_const, 'double');
fclose(fileID);

fprintf("Rho written to files\n")