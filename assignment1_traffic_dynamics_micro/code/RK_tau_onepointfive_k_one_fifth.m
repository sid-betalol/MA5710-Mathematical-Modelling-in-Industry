% Reference used: http://matlab.imm.uran.ru/mirrors/www.cs.runet.edu/~thompson/webddes/tutorial.html
% Reference used: https://youtu.be/TCWrD3cZG9s
% Reference used: https://in.mathworks.com/help/matlab/ref/dde23.html
% Reference used: https://in.mathworks.com/matlabcentral/fileexchange/3899-tutorial-on-solving-ddes-with-dde23
% RK Method
% L = 6m = 0.006 km
% rho_max = 40 cars/km
% d = 19m = 0.019 km
% N = #(cars) = 5
% ts = 1
% vk(ts^2)e = 15.1m = 0.0151 km
% k = 0.2
% tau = 1.5

clear
clc 

global tau 
tau = 1.5; % reaction time of driver
r = 2;

t_span = linspace(0, 40, 1000); % duration of observation
sol = dde23(@func, tau, @zhistory, t_span); %solving the DDe using dde23

% loop for plotting z2,z3,z4,z5
for i=1:4
    d = 0.019*ones(size(sol.y(i,:)));
    y1 = sol.y(i,:) - (i*d);
    plot(sol.x, y1)
    hold on
end

title('Reaction Time = 1.5 seconds and k = 0.2')
xlabel('Time (in seconds)')
ylabel('Shifted perturbation displacements [Zi(t)-(i-1)d] (in kilometers)')
hold off

function Z = func(t, z, zl)
    global tau
    r = 2; % r = 1 if k = 0.1 and r = 2 if k =  0.2 
    Z = [(1/36)*log(1+(40/exp(1))*(-(((15.1*r)/2000)*(1-(t-tau+1)*exp(-t+tau)))-zl(1))) % this row represents z2 
         (1/36)*log(1+(40/exp(1))*(zl(1)-zl(2))) % this row represents z3
         (1/36)*log(1+(40/exp(1))*(zl(2)-zl(3))) % this row represents z4 
         (1/36)*log(1+(40/exp(1))*(zl(3)-zl(4))) % this row represents z5 
        ];
end

% history function
function z = zhistory(t)
    z = [0 0 0 0];
end