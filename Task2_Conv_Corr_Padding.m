%% EE3005 Image Processing - Assignment 2
%  Task 2: Convolution, correlation and border effects
clear; close all; clc;
thisDir = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'helpers'));
fprintf('\n================ TASK 2: CONVOLUTION / CORRELATION / PADDING ================\n');

%% 2.1 Kernel operations on the 5x5 matrix of Figure 1(a)
A = [10 20 30 40 50;
     15 25 35 45 55;
     20 30 40 50 60;
     25 35 45 55 65;
     30 40 50 60 70];
K = [1 0 -1;
     2 0 -2;
     1 0 -1];

Nb = A(2:4, 2:4);                 % 3x3 neighbourhood of the centre pixel A(3,3)
fprintf('Neighbourhood of the centre pixel A(3,3) = %d:\n', A(3,3)); disp(Nb);

% correlation: kernel used as it is
Pcorr = Nb .* K;
corrVal = sum(Pcorr(:));
fprintf('Correlation - element-wise products N.*K:\n'); disp(Pcorr);
fprintf('Correlation result at centre = %d\n\n', corrVal);

% convolution: kernel flipped horizontally and vertically (180 deg rotation)
Kflip = rot90(K, 2);              % same as flipud(fliplr(K))
fprintf('Flipped kernel rot90(K,2):\n'); disp(Kflip);
Pconv = Nb .* Kflip;
convVal = sum(Pconv(:));
fprintf('Convolution - element-wise products N.*Kflip:\n'); disp(Pconv);
fprintf('Convolution result at centre = %d\n', convVal);
fprintf('Check: rot90(K,2) == -K ? %d (K is anti-symmetric left-right)\n\n', isequal(Kflip, -K));

% check with built-in functions (zero padding)
Rcorr = imfilter(A, K, 'corr', 0);
Rconv = imfilter(A, K, 'conv', 0);
Rconv2 = conv2(A, K, 'same');
fprintf('imfilter(A,K,''corr'')  full 5x5 output (zero padding):\n'); disp(Rcorr);
fprintf('imfilter(A,K,''conv'')  full 5x5 output (zero padding):\n'); disp(Rconv);
fprintf('Centre values: imfilter corr = %d, imfilter conv = %d, conv2 = %d\n', ...
        Rcorr(3,3), Rconv(3,3), Rconv2(3,3));

%% 2.2 Border padding with a 3x3 arithmetic mean filter
I = readGray('chest_xray.png');
h = fspecial('average', 3);

out = cell(1,4);
out{1} = imfilter(I, h, 0);             % zero padding
out{2} = imfilter(I, h, 'replicate');   % replication (clamp)
out{3} = imfilter(I, h, 'symmetric');   % reflection (mirror)
out{4} = imfilter(I, h, 'circular');    % circular (wrap)
names = {'Zero', 'Replicate (clamp)', 'Symmetric (mirror)', 'Circular (wrap)'};

fig = newFig('Task 2.2 - padding, full images', [40 40 1400 420]);
for k = 1:4
    subplot(1,4,k); imshow(out{k}); title(['3x3 mean, ' names{k}]);
end
saveFig(fig, 'Fig2_1_padding_full');

% zoom on two corners: top-left (12x12) and bottom-right (12x12)
[M, N] = size(I);
cTL = {1:12, 1:12};  cBR = {M-11:M, N-11:N};
fig = newFig('Task 2.2 - padding, corners', [40 40 1400 700]);
subplot(2,5,1); imshow(I(cTL{:}), [0 1]); title('Original top-left');
subplot(2,5,6); imshow(I(cBR{:}), [0 1]); title('Original bottom-right');
for k = 1:4
    subplot(2,5,k+1); imshow(out{k}(cTL{:}), [0 1]); title(names{k});
    subplot(2,5,k+6); imshow(out{k}(cBR{:}), [0 1]); title(names{k});
end
saveFig(fig, 'Fig2_2_padding_corners');

% difference to the original along the border frame (1 pixel wide)
frame = false(M, N); frame([1 M], :) = true; frame(:, [1 N]) = true;
inner = false(M, N); inner(2:M-1, 2:N-1) = true;
fprintf('\n2.2 Border analysis, 3x3 mean filter (intensity 0-1)\n');
fprintf('Original: mean of first row = %.3f, last row = %.3f, first col = %.3f, last col = %.3f\n', ...
        mean(I(1,:)), mean(I(M,:)), mean(I(:,1)), mean(I(:,N)));
fprintf('%-20s %8s %8s %8s %8s %12s %12s\n', 'Padding', 'row 1', 'row M', 'col 1', 'col N', ...
        'MAE frame', 'MAE inner');
for k = 1:4
    d = abs(out{k} - I);
    fprintf('%-20s %8.3f %8.3f %8.3f %8.3f %12.4f %12.4f\n', names{k}, mean(out{k}(1,:)), ...
            mean(out{k}(M,:)), mean(out{k}(:,1)), mean(out{k}(:,N)), mean(d(frame)), mean(d(inner)));
end
fprintf('(MAE = mean absolute difference to the original image)\n');
fprintf('Corner pixel (1,1): original %.3f, zero %.3f, replicate %.3f, symmetric %.3f, circular %.3f\n', ...
        I(1,1), out{1}(1,1), out{2}(1,1), out{3}(1,1), out{4}(1,1));
fprintf('Corner pixel (M,N): original %.3f, zero %.3f, replicate %.3f, symmetric %.3f, circular %.3f\n', ...
        I(M,N), out{1}(M,N), out{2}(M,N), out{3}(M,N), out{4}(M,N));

% extra: a large 31x31 mean filter makes the border effect clearly visible
h31 = fspecial('average', 31);
big = {imfilter(I, h31, 0), imfilter(I, h31, 'replicate'), ...
       imfilter(I, h31, 'symmetric'), imfilter(I, h31, 'circular')};
fig = newFig('Task 2.2 - padding, 31x31 mean', [40 40 1400 420]);
for k = 1:4
    subplot(1,4,k); imshow(big{k}); title(['31x31 mean, ' names{k}]);
end
saveFig(fig, 'Fig2_3_padding_31x31');
fprintf('\n31x31 mean: mean of the outer 15-pixel band (original = %.3f)\n', ...
        mean(I(~padarray(true(M-30, N-30), [15 15], false))));
band = ~padarray(true(M-30, N-30), [15 15], false);
for k = 1:4
    fprintf('   %-20s %.3f\n', names{k}, mean(big{k}(band)));
end
