clear
clc 

% Initialising variables
x = linspace(-4, 4, 81);   % input range of x values in steps of 0.1 
t = linspace(0, 3, 121);   % input range of t values in steps of 0.025
k = t(2) - t(1);    % Time step size
h = x(2) - x(1);    % Space step size
Nt = 121; % number of time steps
Nx = 81;  % number of space steps
j = @(ro) 8 * ro - 4 * ro.^2; % defining flux j(rho)
j_dash = @(ro) 8 - 8 * ro;    % defining j'(rho)
RHO = zeros(length(x), length(t));

% Initial conditions
RHO(:, 1) = 1 * (x <= 1) + 0.5 * ((x > 1) & (x <= 3)) + 1.5 * (x > 3);
% Boundary conditions
RHO(1, :) = 1;
RHO(81, :) = 1.5;

% Implementing the Gudonov numerical scheme
for idt = 1 : Nt - 1
    for idx = 2 : Nx - 1
        % For computing (rho_i)*
        if (j_dash(RHO(idx, idt)) >= 0) & (j_dash(RHO(idx + 1, idt)) >= 0)
            RHO_star2 = RHO(idx, idt);            
        elseif (j_dash(RHO(idx, idt)) < 0) & (j_dash(RHO(idx + 1, idt)) < 0) 
            RHO_star2 = RHO(idx + 1, idt);            
        elseif (j_dash(RHO(idx, idt)) >= 0) & (j_dash(RHO(idx + 1, idt)) < 0)
            si1 = (j(RHO(idx + 2, idt)) - j(RHO(idx, idt))) / (RHO(idx + 2, idt) - RHO(idx, idt));
            if si1 >= 0
                RHO_star2 = RHO(idx, idt);
            else
                RHO_star2 = RHO(idx + 1, idt);
            end            
        elseif (j_dash(RHO(idx, idt)) < 0) & (j_dash(RHO(idx + 1, idt)) >= 0)
            syms rho
            eqn = j_dash(rho) == 0;
            RHO_star2 = solve(eqn, rho);
        end
        
        % For computing (rho_i-1)*
        if (j_dash(RHO(idx -1, idt)) >= 0) & (j_dash(RHO(idx, idt)) >= 0)
            RHO_star1 = RHO(idx - 1, idt);            
        elseif (j_dash(RHO(idx - 1, idt)) < 0) & (j_dash(RHO(idx, idt)) < 0)
            RHO_star1 = RHO(idx, idt);            
        elseif (j_dash(RHO(idx - 1, idt)) >= 0) & (j_dash(RHO(idx, idt)) < 0)
            si = (j(RHO(idx + 1, idt)) - j(RHO(idx - 1, idt))) / (RHO(idx + 1, idt) - RHO(idx - 1, idt));
            if si >= 0
                RHO_star1 = RHO(idx - 1, idt);
            else
                RHO_star1 = RHO(idx, idt);
            end
        elseif (j_dash(RHO(idx - 1, idt)) < 0) & (j_dash(RHO(idx, idt)) >= 0)
            syms rho
            eqn = j_dash(rho) == 0;     % If j is changed, this can still be used
            RHO_star1 = solve(eqn, rho);
        end
        
        RHO(idx, idt + 1) = RHO(idx, idt) - (k / h) * (j(RHO_star2) - j(RHO_star1));
    end
end

% Plotting the figure
figure(1)
surf(t,x,RHO)
xlabel('$t$','interpreter','latex', 'FontSize', 15)
ylabel('$x$','interpreter','latex', 'FontSize', 15)
zlabel('$\rho(x,t)$','interpreter','latex', 'FontSize', 15)
title('Green Light Problem: Gudonov Scheme', 'FontSize', 12)
colorbar
set(gca, 'XDir','reverse')