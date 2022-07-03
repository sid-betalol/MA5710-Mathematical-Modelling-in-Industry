clear
clc 

% Initialising variables
x = linspace(-4, 4, 81); % input range of x values in steps of 0.1
t = linspace(0, 3, 61);  % input range of t values in steps of 0.1
k = t(2) - t(1);    % Time step size
h = x(2) - x(1);    % Space step size
Nt = 61;    % number of time steps
Nx = 81;    % number of space steps
j = @(ro) ro - ro.^2;   % defining flux j(rho)
j_dash = @(ro) 1 - 2 * ro;  % defining j'(rho)
rho = zeros(length(x), length(t));

% Initial conditions
rho(:, 1) = 1 * (x <= 0) + 0 * (x > 0);  % Density at t=0
% Boundary conditions
rho(1, :) = 1;
rho(81, :) = 0;

% Implementing the Upwind numerical scheme
for idt = 1 : Nt - 1    
    for idx = 2 : Nx - 1
        if j_dash(rho(idx, idt)) > 0
            rho(idx, idt + 1) = rho(idx, idt) - (k / h) * (j(rho(idx, idt)) - j(rho(idx - 1, idt)));  
        elseif j_dash(rho(idx, idt)) < 0
            rho(idx, idt + 1) = rho(idx, idt) - (k / h) * (j(rho(idx + 1, idt)) - j(rho(idx, idt)));
        end    
    end
end

% Plotting the figure
figure(1)
surf(t,x,rho)
xlabel('$t$','interpreter','latex', 'FontSize', 15)
ylabel('$x$','interpreter','latex', 'FontSize', 15)
zlabel('$\rho(x,t)$','interpreter','latex', 'FontSize', 15)
colorbar
title('Green Light Problem: Upwind Scheme', 'FontSize', 12)
set(gca, 'XDir','reverse')