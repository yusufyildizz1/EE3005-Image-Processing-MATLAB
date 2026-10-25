%% EE3005 Image Processing - Assignment 2
%  Task 7: Edge detection with Sobel, Laplacian, LoG and DoG
%  Image: axial T1 brain MRI slice (images/brain_mri.png)
clear; close all; clc;
thisDir = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'helpers'));
fprintf('\n================ TASK 7: EDGE DETECTION ================\n');

I = readGray('brain_mri.png');
rng(7); In = min(max(I + 0.02*randn(size(I)), 0), 1);   % noisy copy, sigma = 0.02
lapK = [0 1 0; 1 -4 1; 0 1 0];                          % Laplacian kernel (lecture)

%% helper: all responses for one input image
respond = @(X) struct( ...
    'sobel',  sobelMag(X), ...
    'lap',    imfilter(X, lapK, 'replicate'), ...
    'log1',   imfilter(imfilter(X, gaussKernel(1), 'replicate'), lapK, 'replicate'), ...
    'log2',   imfilter(imfilter(X, gaussKernel(2), 'replicate'), lapK, 'replicate'), ...
    'dog1',   imfilter(X, gaussKernel(1), 'replicate') - imfilter(X, gaussKernel(1.6), 'replicate'), ...
    'dog2',   imfilter(X, gaussKernel(2), 'replicate') - imfilter(X, gaussKernel(3.2), 'replicate'), ...
    'dog3',   imfilter(X, gaussKernel(1), 'replicate') - imfilter(X, gaussKernel(4), 'replicate'));
Rc = respond(I);     % clean image
Rn = respond(In);    % noisy image

%% 7.1 Sobel and Laplacian
[~, Gx, Gy] = sobelMag(I);
bwSobel = edge(I, 'sobel');
fig = newFig('Task 7.1 - Sobel and Laplacian', [40 40 1400 750]);
subplot(2,4,1); imshow(I);               title('Brain MRI (T1, axial)');
subplot(2,4,2); imshow(Gx, []);          title('Sobel G_x (vertical edges)');
subplot(2,4,3); imshow(Gy, []);          title('Sobel G_y (horizontal edges)');
subplot(2,4,4); imshow(Rc.sobel, []);    title('Sobel magnitude');
subplot(2,4,5); imshow(bwSobel);         title('edge(I,''sobel'') binary');
subplot(2,4,6); imshow(Rc.lap, []);      title('Laplacian response (signed)');
subplot(2,4,7); imshow(abs(Rc.lap), []); title('|Laplacian|');
subplot(2,4,8); imshow(Rn.lap, []);      title('Laplacian on noisy image');
saveFig(fig, 'Fig7_1_sobel_laplacian');

% 1-D profile (head border + lateral ventricles) to compare 1st and 2nd derivative
pr = 85; pc = 10:120;                     % crosses the head border and the lateral ventricles
fig = newFig('Task 7.1 - profile', [40 40 1100 650]);
subplot(3,1,1); plot(pc, I(pr,pc), 'k', 'LineWidth', 1.5); grid on; ylabel('intensity');
title(sprintf('Row %d: intensity, first and second derivative', pr));
subplot(3,1,2); plot(pc, Gx(pr,pc), 'b', 'LineWidth', 1.2); grid on; ylabel('Sobel G_x');
subplot(3,1,3); plot(pc, Rc.lap(pr,pc), 'r', 'LineWidth', 1.2); hold on;
plot(pc, zeros(size(pc)), 'k:'); hold off; grid on; ylabel('Laplacian'); xlabel('column');
saveFig(fig, 'Fig7_2_derivative_profile');

%% 7.2 LoG
bwLog1 = edge(I, 'log', [], 1);
bwLog2 = edge(I, 'log', [], 2);
fig = newFig('Task 7.2 - LoG', [40 40 1400 750]);
subplot(2,3,1); imshow(Rc.log1, []); title('LoG response, \sigma = 1');
subplot(2,3,2); imshow(Rc.log2, []); title('LoG response, \sigma = 2');
subplot(2,3,3); imshow(Rn.log2, []); title('LoG \sigma = 2 on noisy image');
subplot(2,3,4); imshow(bwLog1);      title('zero crossings, \sigma = 1');
subplot(2,3,5); imshow(bwLog2);      title('zero crossings, \sigma = 2');
subplot(2,3,6); imshow(edge(In, 'log', [], 2)); title('zero crossings \sigma = 2, noisy');
saveFig(fig, 'Fig7_3_LoG');

%% 7.3 DoG
fig = newFig('Task 7.3 - DoG', [40 40 1400 750]);
subplot(2,3,1); imshow(Rc.dog1, []); title('DoG \sigma_1=1, \sigma_2=1.6');
subplot(2,3,2); imshow(Rc.dog2, []); title('DoG \sigma_1=2, \sigma_2=3.2');
subplot(2,3,3); imshow(Rc.dog3, []); title('DoG \sigma_1=1, \sigma_2=4');
subplot(2,3,4); imshow(Rn.dog1, []); title('noisy: \sigma_1=1, \sigma_2=1.6');
subplot(2,3,5); imshow(Rn.dog2, []); title('noisy: \sigma_1=2, \sigma_2=3.2');
subplot(2,3,6); imshow(Rn.dog3, []); title('noisy: \sigma_1=1, \sigma_2=4');
saveFig(fig, 'Fig7_4_DoG');

% DoG <-> LoG relationship: correlation with LoG at the middle scale
fprintf('7.3 DoG vs LoG (correlation coefficient over the image)\n');
pairs = [1 1.6; 2 3.2; 1 4];
for k = 1:3
    s1 = pairs(k,1); s2 = pairs(k,2);
    dog = imfilter(I, gaussKernel(s1), 'replicate') - imfilter(I, gaussKernel(s2), 'replicate');
    sm  = sqrt(s1*s2);
    lg  = imfilter(imfilter(I, gaussKernel(sm), 'replicate'), lapK, 'replicate');
    fprintf('   DoG(%.1f, %.1f) vs LoG(sigma = %.2f): r = %+.3f\n', s1, s2, sm, corr2(dog, lg));
end
fprintf('   (negative sign: G_s1 - G_s2 with s1 < s2 is the inverted LoG)\n');

%% 7.4 Comparison: noise sensitivity and cost
% Signal: strongest 1 % of the |response| on the CLEAN image (real structures).
% Noise : std of the response in the background of the NOISY image (outside
%         the head the clean image is 0, so any response there is noise).
bg = imerode(I < 0.02, ones(11));           % background, >= 5 px away from the head
fprintf('\n7.4 Edge-to-noise ratio, ENR = P99(|response clean|) / RMS(response noisy, background)\n');
names = {'sobel', 'lap', 'log1', 'log2', 'dog1', 'dog2', 'dog3'};
label = {'Sobel magnitude', 'Laplacian', 'LoG s=1', 'LoG s=2', 'DoG 1/1.6', 'DoG 2/3.2', 'DoG 1/4'};
mult  = {'2x9 = 18', '9', '49+9 = 58', '169+9 = 178', '49+121 = 170', '169+441 = 610', '49+625 = 674'};
multS = {'2x6 = 12', '9', '14+9 = 23', '26+9 = 35', '14+22 = 36', '26+42 = 68', '14+50 = 64'};
fprintf('%-16s %8s %16s %18s\n', 'Method', 'ENR', 'mult/pixel 2-D', 'mult/pixel separ.');
for k = 1:numel(names)
    rn = Rn.(names{k}); rc = abs(Rc.(names{k}));
    sv = sort(rc(:)); sig = sv(round(0.99*numel(sv)));
    enr = sig / sqrt(mean(rn(bg).^2));          % RMS noise response in background
    fprintf('%-16s %8.1f %16s %18s\n', label{k}, enr, mult{k}, multS{k});
end
fprintf('(Gaussian kernels are (2*ceil(3*sigma)+1)^2; separable = two 1-D passes)\n');
% timing of each operator alone
fprintf('\nTiming per operator (mean of 20 runs):\n');
ops = {@() sobelMag(I), @() imfilter(I, lapK, 'replicate'), ...
       @() imfilter(imfilter(I, gaussKernel(2), 'replicate'), lapK, 'replicate'), ...
       @() imfilter(I, gaussKernel(2), 'replicate') - imfilter(I, gaussKernel(3.2), 'replicate')};
opN = {'Sobel', 'Laplacian', 'LoG s=2', 'DoG 2/3.2'};
for k = 1:4
    ops{k}(); tic; for r = 1:20, ops{k}(); end
    fprintf('   %-10s %.2f ms\n', opN{k}, 1000*toc/20);
end

fig = newFig('Task 7.4 - comparison on noisy image', [40 40 1400 420]);
cm = {Rn.sobel, abs(Rn.lap), abs(Rn.log2), abs(Rn.dog2)};
cn = {'Sobel |grad|', '|Laplacian|', '|LoG| \sigma=2', '|DoG| 2/3.2'};
for k = 1:4
    subplot(1,4,k); imshow(cm{k}, []); title([cn{k} ' (noisy)']);
end
saveFig(fig, 'Fig7_5_comparison_noisy');
