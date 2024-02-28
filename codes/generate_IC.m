close all; clc; clear;
%% Setting up random seed
rng(2022);

%% Setting up simulation parameters
Params

%% Initializing order parameter n
n = zeros([Nx,Ny],'quaternion');
rand_pos = zeros(np, 2);
qlist = zeros(np,1,'quaternion');
nearest = zeros(Nx,Ny);

curr_grain_id = 1;
rand_pos(:,1) = ceil(rand(np,1).*Nx);
rand_pos(:,2) = ceil(rand(np,1).*Ny);
while curr_grain_id <= np
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

%% Write n to files
save('../data/initial_n.mat',"n");