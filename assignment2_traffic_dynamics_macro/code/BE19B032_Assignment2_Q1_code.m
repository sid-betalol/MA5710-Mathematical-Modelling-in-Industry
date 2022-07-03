clear
clc

x = linspace(-4, 4, 81);    % input range of x values in steps of 0.1
t = linspace(0, 3, 31);     % input range of t values in steps of 0.1     
[X,T] = meshgrid(x,t);
rho = zeros(length(t), length(x)); 

% different cases
case1 = X <= -T;
case2 = (-T < X) & (X <= T);
case3 = X > T;

% density value corresponding to each case
rho(case1) = 1;
rho(case2) = (0.5) * (1 - (X(case2) ./ T(case2)));
rho(case3) = 0;

% plotting the figure
figure(1)
surf(x, t, rho)
ylabel('$t$','interpreter','latex', 'FontSize', 15)
xlabel('$x$','interpreter','latex', 'FontSize', 15)
zlabel('$\rho(x,t)$','interpreter','latex', 'FontSize', 15)
title('Green Light Problem: Exact Solution',  'FontSize', 15)
colorbar