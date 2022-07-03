% references used: Prof. Sundar's classnotes
% references used: https://staff.fnwi.uva.nl/r.vandenboomgaard/nldiffusionweb/nldiffusioncode.pdf

clear
clc

% storing and displaying original image
img_1 = imread('cameraman.tif');
imshow(img_1);

% storing and displaying noisy version of above image
img_2 = im2double(imnoise(img_1,'gaussian', 0, 0.01));
figure(1)
imshow(img_2)

% Perona-Malik filtering for various contrast parameter values
l = [1 0.5 0.1 0.01 0.005 0.001 0.0001];
for i = 1:7
    img_final = PM(img_2, 1.6, 20, l(i));
    disp(['(Lambda: ', num2str(l(i)), ') PSNR=', num2str(psnr(img_final, im2double(img_1)))]);
    figure(1+i)
    imshow(img_final)
end

% estimating the contrast parameter
[x_, y_] = gradient(img_2); % function for finding gradient
G = x_.^2 + y_.^2; % magnitude square
[counts, bins] = imhist(G(:));  % histogram of gradient matrix
cdf = cumsum(counts);   
ncdf = cdf / cdf(size(bins, 1)); 
idx = find(ncdf >= 0.95);
l_best = bins(min(idx))    
img_final = PM(img_2, 1.7, 20, l_best);
figure
imshow(img_final)
disp(['For l_best, PSNR=', num2str(psnr(img_final, im2double(img_1)))])

% finding stopping criteria
for t=0.5:0.1:3
    img_final = PM(img_2, t, 20, l_best);
    disp(['(t=', num2str(t), ') PSNR=' ,num2str(psnr(img_final, im2double(img_1)))])
end 

% function for applying Perona-Malik filter
function img_final = PM(img, t, N_iter, l)
    dt = t / N_iter; % time step size
    U = zeros(size(img,1), size(img,2), N_iter+1); % array for storing our image after each time step
    U(:, :, 1) = img; % initial condition
    C = zeros(size(img,1), size(img,2), N_iter+1);
    
    % finding C at t=0
    C(2:end-1, 2:end-1, 1) = ...
        ones(size(img,1)-2, size(img,2)-2) ./ (ones(size(img,1)-2, size(img,2)-2) ...
        + (0.25/l^2)*((U(3:end, 2:end-1, 1) - U(1:end-2, 2:end-1, 1)).^2 ...
        + (U(2:end-1, 3:end, 1) - U(2:end-1, 1:end-2, 1)).^2));
    C(1, 2:end-1, 1) = ...
        ones(1, size(img,2)-2) ./ (ones(1, size(img,2)-2) + (1/l^2)*((U(2, 2:end-1, 1) - U(1, 2:end-1, 1)).^2 ...
        + (U(1, 3:end, 1) - U(1, 2:end-1, 1)).^2));      
    C(end, 2:end-1, 1) = ones(1, size(img,2)-2) ./ (ones(1, size(img,2)-2) + (1/l^2)*((U(end, 2:end-1, 1) - U(end-1, 2:end-1, 1)).^2 ...
        + (U(end, 3:end, 1) - U(end, 2:end-1, 1)).^2));  
    C(2:end-1, 1, 1) = ...
        ones(size(img, 1)-2, 1) ./ (ones(size(img, 1)-2, 1) + (1/l^2)*((U(3:end, 1, 1) - U(2:end-1, 1, 1)).^2 ...
        + (U(2:end-1, 2, 1) - U(2:end-1, 1, 1)).^2));    
    C(2:end-1, end, 1) = ...
        ones(size(img, 1)-2, 1) ./ (ones(size(img, 1)-2, 1) + (1/l^2)*((U(3:end, end, 1) - U(2:end-1, end, 1)).^2 ...
        + (U(2:end-1, end-1, 1) - U(2:end-1, end, 1)).^2));  

    for i = 1:N_iter
        
        % calculate image at next time step
        U(2:end-1, 2:end-1, i+1) = ...
            U(2:end-1, 2:end-1, i) + 0.5*dt*...
            ((C(3:end, 2:end-1, i) + C(2:end-1, 2:end-1, i)).*(U(3:end, 2:end-1, i) - U(2:end-1, 2:end-1, i)) -...
            (C(2:end-1, 2:end-1, i) + C(1:end-2, 2:end-1, i)).*(U(2:end-1, 2:end-1, i) - U(1:end-2, 2:end-1, i)) +...
            (C(2:end-1, 3:end, i) + C(2:end-1, 2:end-1, i)).*(U(2:end-1, 3:end, i) - U(2:end-1, 2:end-1, i)) -...
            (C(2:end-1, 2:end-1, i) + C(2:end-1, 1:end-2, i)).*(U(2:end-1, 2:end-1, i) - U(2:end-1, 1:end-2, i)));
        
        % enforce boundary conditions
        U(1, 2:end-1, i+1) = U(2, 2:end-1, i+1);
        U(end, 2:end-1, i+1) = U(end-1, 2:end-1, i+1);
        U(2:end-1, 1, i+1) = U(2:end-1, 2, i+1);
        U(2:end-1, end, i+1) = U(2:end-1, end-1, i+1);
        U(1, 1, i+1) = 0;
        U(1, end, i+1) = 0;
        U(end, 1, i+1) = 0;
        U(end, end, i+1) = 0;
        
        % finding C at next time step
        C(2:end-1, 2:end-1, i+1) = ...
            ones(size(img,1)-2, size(img,2)-2) ./ (ones(size(img,1)-2, size(img,2)-2) ...
            + (0.25/l^2)*((U(3:end, 2:end-1, i+1) - U(1:end-2, 2:end-1, i+1)).^2 ...
            + (U(2:end-1, 3:end, i+1) - U(2:end-1, 1:end-2, i+1)).^2));
        C(1, 2:end-1, i+1) = ...
            ones(1, size(img,2)-2) ./ (ones(1, size(img,2)-2) + ...
            (1/l^2)*(U(1, 3:end, i+1) - U(1, 2:end-1, i+1)).^2);
        C(end, 2:end-1, i+1) = ...
            ones(1, size(img,2)-2) ./ (ones(1, size(img,2)-2) + ...
            (1/l^2)*(U(end, 3:end, i+1) - U(end, 2:end-1, i+1)).^2);
        C(2:end-1, 1, i+1) = ...
            ones(size(img, 1)-2, 1) ./ (ones(size(img, 1)-2, 1) + ...
            (1/l^2)*(U(3:end, 1, i+1) - U(2:end-1, 1, i+1)).^2);
        C(2:end-1, end, i+1) = ...
            ones(size(img, 1)-2, 1) ./ (ones(size(img, 1)-2, 1) + ...
            (1/l^2)*(U(3:end, end, i+1) - U(2:end-1, end, i+1)).^2);
    end
    
    img_final = U(:, :, end); % final image at time t
end