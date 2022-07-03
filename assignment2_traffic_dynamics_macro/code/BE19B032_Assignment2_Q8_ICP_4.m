clear
clc 

% Initialising variables
x = linspace(-4, 4, 51);    % input range of x values in steps of 0.16
t = linspace(0, 3, 121);    % input range of t values in steps of 0.025
k = t(2) - t(1);    % Time step size
h = x(2) - x(1);    % Space step size
Nt = 121; % number of time steps
Nx = 51;  % number of space steps
j = @(ro) 8 * ro - 4 * ro.^2; % defining flux j(rho)
j_dash = @(ro) 8 - 8 * ro; % defining j'(rho)
rho = zeros(length(x), length(t));

% Initial conditions
rho(:, 1) = 1 * (x <= 1) + 0.5 * ((x > 1) & (x <= 3)) + 1.5 * (x > 3);
% Boundary conditions
rho(1, :) = 1;
rho(51, :) = 1.5;

% Implementing the Mac-Cormack numerical scheme
for idt = 1 : Nt - 1
    for idx = 2 : Nx - 1
        rho_star1 = rho(idx - 1, idt) - (k / h) * (j(rho(idx, idt)) - j(rho(idx - 1, idt)));
        rho_star2 = rho(idx, idt) - (k / h) * (j(rho(idx + 1, idt)) - j(rho(idx, idt)));
        rho(idx, idt + 1) = 0.5 * (rho(idx, idt) + rho_star2) - (k / h) * (j(rho_star2) - j(rho_star1));
    end
end

% Plotting the figure
figure(1)
surf(t,x,rho)
xlabel('$t$','interpreter','latex', 'FontSize', 15)
ylabel('$x$','interpreter','latex', 'FontSize', 15)
zlabel('$\rho(x,t)$','interpreter','latex', 'FontSize', 15)
title('Intersecting Characteristics: Mac-Cormack Scheme', 'FontSize', 12)
colorbar
set(gca, 'XDir','reverse')