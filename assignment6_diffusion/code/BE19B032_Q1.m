% references used: Prof. Sundar's classnotes
% references used: THE MATHEMATICS. OF DIFFUSION. BY. J. CRANK
% references used: https://in.mathworks.com/matlabcentral/answers/137776-comparing-two-matrix-in-matlab
% the cameraman image from the Image Processing toolbox has been used for all three questions.

clear
clc

% storing and displaying original image
img_1 = imread('cameraman.tif');
imshow(img_1);

% storing and displaying noisy version of above image
img_2 = im2double(imnoise(img_1,'gaussian', 0, 0.01));
figure(1)
imshow(img_2)

% image reconstruction for a list of sigma values along with PSNR
for sigma = 0.1:0.1:1.5
    img_final = linear(img_2, sigma);
    disp(['(sigma=', num2str(sigma), ') PSNR=', num2str(psnr(img_final, im2double(img_1)))])
    figure
    imshow(img_final)
end

% verifying the properties of linear diffusion
% the two images are compared within some tolerance range using the 
% function isequalRel. The tolerance is given by the third argument/input 
% of isequalRel

sigma = 1;
% Translation invariance
disp('Translation invariance-')
A = linear(circshift(img_2, 50, 2), sigma);
B = circshift(linear(img_2, sigma), 50, 2);
isequalRel(A, B, 0.1);

% Gray level shift
disp('Gray level shift-')
A = linear(img_2 + ones(size(img_2)), sigma);
B = linear(img_2, sigma) + ones(size(img_2));
isequalRel(A, B, 0.05);

% Scale invariance
disp('Scale invariance-')
A = linear(imresize(img_2, 5), sigma);
B = imresize(linear(img_2, 5*sigma), 5);
isequalRel(A, B, 0.1);

% Conservation of average value
disp('Conservation of average value-')
A = linear(mean(img_2), sigma);
B = mean(linear(img_2, sigma));
isequalRel(A, B, 0.05);

% Semigroup property
disp('Semigroup property:')
t = 1;
s = 2;
A = linear(img_2, sqrt(2*(s+t)));
x = linear(img_2, sqrt(2*s));
B = linear(x, sqrt(2*t));
isequalRel(A, B, 0.1);

% Isometry invariance
disp('Isometry invariance-')
A = linear(img_2.', sigma);
B = linear(img_2, sigma).';
isequalRel(A, B, 0.05);

% Comparison Principle
disp('Comparison Principle-')
A = linear(img_2, sigma);
B = linear(img_2 + 5*ones(size(img_2)), sigma);
x = B > A;
if isequal(x, ones(size(x)))
    disp('True')
else
    disp('False')
end 

function img_final = linear(img, sigma)
    img_final = imgaussfilt(img, sigma); 
end

function test = isequalRel(x,y,tol)
    test = ( abs(x-y) <= ( tol*max(abs(x),abs(y)) + eps) );
    if all(test)
        disp("True")
    else
        disp("False");
    end
end
