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

% Catte et al. filtering for various contrast parameter values 
sigma = 1; % fix sigma
l = [1 0.5 0.1 0.01 0.005 0.001 0.0001];
for i = 1:7
    img_final = catte(img_2, 1.7, 20, l(i), sigma);
    disp(['(Lambda: ', num2str(l(i)), ') PSNR=', num2str(psnr(img_final, im2double(img_1)))]);
    figure(1+i)
    imshow(img_final)
end

% finding stopping criteria
for t=0.5:0.1:3
    img_final = catte(img_2, t, 20, 0.08, sigma);
    disp(['(t=', num2str(t), ') PSNR=' ,num2str(psnr(img_final, im2double(img_1)))])
end 

% function for applying Catte et al. filter
function img_final = catte(img, t, N_iter, l, sigma)
    dt = t / N_iter;
    U = zeros(size(img,1), size(img,2), N_iter+1);
    U(:, :, 1) = img;
    C = zeros(size(img,1), size(img,2), N_iter+1);

    for i = 1:N_iter
        
        % finding C at current time step
        U_sig = imgaussfilt(U(:, :, i), sigma);
        C(2:end-1, 2:end-1, i) = ...
            ones(size(img,1)-2, size(img,2)-2) ./ (ones(size(img,1)-2, size(img,2)-2) ...
            + (0.25/l^2)*((U_sig(3:end, 2:end-1) - U_sig(1:end-2, 2:end-1)).^2 ...
            + (U_sig(2:end-1, 3:end) - U_sig(2:end-1, 1:end-2)).^2));
        C(1, 2:end-1, i) = ...
            ones(1, size(img,2)-2) ./ (ones(1, size(img,2)-2) + ...
            (1/l^2)*((U_sig(2, 2:end-1) - U_sig(1, 2:end-1)).^2 ...
            + (U_sig(1, 3:end) - U_sig(1, 2:end-1)).^2));      
        C(end, 2:end-1, i) = ...
            ones(1, size(img,2)-2) ./ (ones(1, size(img,2)-2) + ...
            (1/l^2)*((U_sig(end, 2:end-1) - U_sig(end-1, 2:end-1)).^2 ...
            + (U_sig(end, 3:end) - U_sig(end, 2:end-1)).^2)); 
        C(2:end-1, 1, i) = ...
            ones(size(img, 1)-2, 1) ./ (ones(size(img, 1)-2, 1) + ...
            (1/l^2)*((U_sig(3:end, 1) - U_sig(2:end-1, 1)).^2 ...
            + (U_sig(2:end-1, 2) - U_sig(2:end-1, 1, 1)).^2));    
        C(2:end-1, end, i) = ...
            ones(size(img, 1)-2, 1) ./ (ones(size(img, 1)-2, 1) + ...
            (1/l^2)*((U_sig(3:end, end) - U_sig(2:end-1, end)).^2 ...
            + (U_sig(2:end-1, end-1) - U_sig(2:end-1, end)).^2)); 
    
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
    end
    
    img_final = U(:, :, end); % final image at time t
end