dx = 0.01;
x0 = 0;
lx = 2.00;

xf = lx - dx;

np = 2;

Nx = round((xf - x0) / dx + 1);

l_gb_ratio = 0.1;

kappa = (1/8)*l_gb_ratio^2;%0.00125;
m0 = 1;
L = 1;
alpha = 1.5;
beta = (2+pi)/(4-pi) / 8;

%% Material parameters
rho_mean = 1.21e13;
rho_std = 3.8e12;
rho_0 = 7.0e12; % base dislocation density for swept grains 

l0_const = 128e-6;
gamma_gb = 0.595;
m0_const = 6.0 * gamma_gb / (l_gb_ratio*l0_const);

W = sqrt(kappa / m0);
G = 28.4e9;
b = 0.255e-9;