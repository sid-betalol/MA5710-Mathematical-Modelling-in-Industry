clear
clc

x = linspace(-4, 4, 81);    % input range of x values in steps of 0.1
t = linspace(0, 3, 31);     % input range of t values in steps of 0.1   
[X,T] = meshgrid(x,t);
rho = zeros(length(t), length(x));

% different cases
case1a = (X <= 1) & (T <= 0.5);
case1b = (X > 1) & (X <= (1 + 4*T)) & (T <= 0.5);
case1c = (X > (1 + 4*T)) & (X <= 3) & (T <= 0.5);
case1d =  (X > 3) & (T <= 0.5);
case2a = (X <= 1) & (T > 0.5) & (T <= 2);
case2b = (X > 1) & (X <= (1 + 4*sqrt(2*T) - 4*T)) & (T > 0.5) & (T <= 2);
case2c = (X > (1 + 4*sqrt(2*T) - 4*T)) & (T > 0.5) & (T <= 2);
case3a = (X <= (5 - 2*T)) & (T > 2);
case3b = (X > (5 - 2*T)) & (T > 2);

% density value corresponding to each case
rho(case1a) = 1;
rho(case1b) = 1 - (1/8)*((X(case1b) ./ T(case1b)) - (1 ./ T(case1b)));
rho(case1c) = 0.5;
rho(case1d) = 1.5;
rho(case2a) = 1;
rho(case2b) = 1 - (1/8)*((X(case2b) ./ T(case2b)) - (1 ./ T(case2b)));
rho(case2c) = 1.5;
rho(case3a) = 1;
rho(case3b) = 1.5;

% plotting the figure
figure(1)
surf(x, t, rho)
ylabel('$t$','interpreter','latex', 'FontSize', 15)
xlabel('$x$','interpreter','latex', 'FontSize', 15)
zlabel('$\rho(x,t)$','interpreter','latex', 'FontSize', 15)
title('Intersecting Characteristics Problem: Exact Solution',  'FontSize', 12)
colorbar