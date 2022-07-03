clear
clc 

% Initialising variables
x = linspace(-4, 4, 41);     % input range of x values in steps of 0.2
t = linspace(0, 3, 151);     % input range of t values in steps of 0.02
k = t(2) - t(1);    % Time step size
h = x(2) - x(1);    % Space step size
Nt = 151;   % number of time steps
Nx = 41;    % number of space steps
j = @(ro) 8 * ro - 4 * ro.^2; % defining flux j(rho)
j_dash = @(ro) 8 - 8 * ro;    % defining j'(rho)
rho = zeros(length(x), length(t));

% Initial conditions
rho(:, 1) = 1 * (x <= 1) + 0.5 * ((x > 1) & (x <= 3)) + 1.5 * (x > 3);
% Boundary conditions
rho(1, :) = 1;
rho(41, :) = 1.5;

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
title('Intersecting Characteristics Problem: Upwind Scheme', 'FontSize', 12)
colorbar
set(gca, 'XDir','reverse')