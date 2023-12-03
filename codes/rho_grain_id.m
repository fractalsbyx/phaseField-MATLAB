close all; clc; clear;
%% Setting up random seed
rng(2021);

%% Params
Params

%% Read in order parameters (grain information)
fileID = fopen('../data/initial_n.dat');
n = fread(fileID, Nx * Ny * np, 'double');
n = reshape(n, [Nx, Ny, np]);
fclose(fileID);

%% Rho Constants
rho_const = normrnd(rho_mean, rho_std, np, 1);

%% Initializaing dislocation density rho
all_sum_nsq = sum(n.^2, 3);
rho = zeros(Nx, Ny);
for i = 1:length(rho_const)
  rho = rho + n(:, :, i).^2 .* rho_const(i);
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