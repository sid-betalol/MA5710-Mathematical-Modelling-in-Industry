% This approach is based on section 8.3.2
% of the book Mathematical Modelling: A case studies approach
% Euler's Method
% k = 0.2
% tau = 0.5

clear
clc

% Initial conditions and setup
K = 200; % 1/K will be the step-size we'll use for Euler's Method  
h = 1;  % step size
n = 8000; % n/K = 40 seconds. The time, we want to observe the model for
M = 100; % M/K = 0.5s , the driver's reaction time
speed = 1/36; % in kilometers per second
t = 0:h:n+M;  % the range of time inputs 
rho = 40; % rho max is given as 40 cars/km.s

% declaring null arrays of Zi's for i from 2 to 5 in which we'll input
% the values of Zi's at differnt values of t
Z2 = zeros(size(t));    
Z3 = zeros(size(t)); 
Z4 = zeros(size(t));  
Z5 = zeros(size(t));  
z1 = @func1;

% For loop to solve the Delay Differential Equation using Euler's method
for i=1:n
    Z2(M+i+1) = (1/K)*(speed)*log(1+(rho/exp(1))*(z1(i-M)-Z2(i))) + Z2(M+i);
end

tf = linspace(0,n+M,n+M+1);
d = (0.019)*ones(size(t)); % d = 19 meters = 0.019 km, given

plot(t,Z2-d)
title('Reaction Time = 0.5 seconds and k = 0.2')
xlabel('time in seconds * 200')
ylabel('shifted perturbation displacements [Zi(t)-(i-1)d] in kilometers')
hold on

for i=1:n
    Z3(M+i+1) = (1/K)*(speed)*log(1+(rho/exp(1))*(Z2(i)-Z3(i))) + Z3(M+i); % y(n+1) = h*f(x,y) + y(n)
end
plot(tf,Z3-2*d)

for i=1:n
    Z4(M+i+1) = (1/K)*(speed)*log(1+(rho/exp(1))*(Z3(i)-Z4(i))) + Z4(M+i);
end
plot(tf,Z4-3*d)

for i=1:n
    Z5(M+i+1) = (1/K)*(speed)*log(1+(rho/exp(1))*(Z4(i)-Z5(i))) + Z5(M+i);
end
plot(tf,Z5-4*d)
hold off

function z = func1(n)
    K = 200;
    if n < 0
        z = 0;
    else 
        % z = -vB(t)
        z = -(0.0151)*(1-((n/K)+1)*exp(-n/K));  %for k = 0.2, vk(ts^2)e = 15.1m = 0.0151 km
    end
end