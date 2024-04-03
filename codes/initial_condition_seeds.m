close all; clc; clear;
%% Setting up random seed
rng(2022);

%% Setting up simulation parameters
Params

%% Generating grid points
center = (xf-x0)/2;
X = x0:dx:xf;

%% Initializing order parameter n
n = zeros(Nx, np);
intf = sqrt(2*kappa/m0);
n(:, 1) = 0.5.*(1-tanh((X-center)./intf));
n(:, 2) = 1 - n(:, 1);

%% Write n to files
fileID = fopen('../data/initial_n.dat','w');
fwrite(fileID, n, 'double');
fclose(fileID);
fprintf("Initial condition written to ../data/initial_n.dat\n")
