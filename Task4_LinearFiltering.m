%% EE3005 Image Processing - Assignment 2
%  Task 4: Linear filtering (arithmetic mean, Gaussian) for noise reduction
%  Reference image: chest X-ray (images/chest_xray.png), double in [0,1]
clear; close all; clc;
thisDir = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'helpers'));
fprintf('\n================ TASK 4: LINEAR FILTERING ================\n');

Iref = readGray('chest_xray.png');        % clean reference, kept unchanged
N = makeNoisyImages(Iref);                % same noisy images are used in Task 5
roiZoom = {150:270, 120:240};             % right lung: vessels and rib edges

%% 4.1 Noise generation
fig = newFig('Task 4.1 - noise', [40 40 1400 750]);
subplot(2,3,1); imshow(Iref);     title('Reference (clean)');
subplot(2,3,2); imshow(N.gauss);  title('Gaussian noise, var = 0.01');
subplot(2,3,3); imshow(N.sp);     title('Salt & pepper, d = 0.05');
subplot(2,3,4); imshow(Iref(roiZoom{:}));    title('Reference - zoom');
subplot(2,3,5); imshow(N.gauss(roiZoom{:})); title('Gaussian - zoom');
subplot(2,3,6); imshow(N.sp(roiZoom{:}));    title('Salt & pepper - zoom');
saveFig(fig, 'Fig4_1_noisy_images');

nSalt = mean(N.sp(:) == 1 & Iref(:) < 1); nPep = mean(N.sp(:) == 0 & Iref(:) > 0);
fprintf('4.1 Gaussian noise: std of (noisy - reference) = %.4f\n', std(N.gauss(:) - Iref(:)));
fprintf('    Salt & pepper: %.2f %% pixels changed to 1, %.2f %% changed to 0\n', 100*nSalt, 100*nPep);

%% 4.2 Mean and Gaussian filters
filt = struct('name', {}, 'h', {});
for k = [3 5 7]
    filt(end+1).name = sprintf('Mean %dx%d', k, k);              %#ok<SAGROW>
    filt(end).h = fspecial('average', k);
end
for sg = [1 2]
    g = gaussKernel(sg);
    filt(end+1).name = sprintf('Gaussian s=%d (%dx%d)', sg, size(g,1), size(g,2)); %#ok<SAGROW>
    filt(end).h = g;
end

noisyList = {N.gauss, N.sp}; noiseName = {'Gaussian', 'Salt & pepper'};
F = cell(numel(filt), 2);
for n = 1:2
    for f = 1:numel(filt)
        F{f,n} = imfilter(noisyList{n}, filt(f).h, 'replicate');
    end
end

for n = 1:2
    fig = newFig(['Task 4.2 - ' noiseName{n}], [40 40 1400 800]);
    subplot(2,3,1); imshow(noisyList{n}); title([noiseName{n} ' noisy']);
    for f = 1:numel(filt)
        subplot(2,3,f+1); imshow(F{f,n}); title(filt(f).name);
    end
    saveFig(fig, sprintf('Fig4_%d_filtered_%s', n+1, lower(strrep(noiseName{n}, ' & ', '_'))));
end

% zoom comparison (Gaussian noise case)
fig = newFig('Task 4.2 - zoom', [40 40 1400 420]);
zImgs = {Iref, N.gauss, F{1,1}, F{3,1}, F{4,1}, F{5,1}};
zName = {'Reference', 'Gaussian noisy', filt(1).name, filt(3).name, filt(4).name, filt(5).name};
for k = 1:6
    subplot(1,6,k); imshow(zImgs{k}(roiZoom{:})); title(zName{k});
end
saveFig(fig, 'Fig4_4_zoom_lung');

%% 4.3 PSNR table
% Blur-only loss: the same filter applied to the CLEAN image. This shows
% how much of the structure the filter removes even without noise.
% Edge retention: mean Sobel magnitude on the strongest 10 % edges of the
% clean image, after / before filtering the clean image.
gRef = sobelMag(Iref);
sv = sort(gRef(:)); edgeMask = gRef >= sv(round(0.9*numel(sv)));

fprintf('\n4.3 PSNR [dB] against the clean reference\n');
fprintf('%-24s %10s %15s %18s %16s\n', 'Image', 'Gaussian', 'Salt & pepper', ...
        'Clean->filter PSNR', 'Edge retention %');
fprintf('%-24s %10.2f %15.2f %18s %16s\n', 'Noisy (no filter)', ...
        computePSNR(Iref, N.gauss), computePSNR(Iref, N.sp), '-', '100.0');
P = zeros(numel(filt), 2);
for f = 1:numel(filt)
    P(f,1) = computePSNR(Iref, F{f,1});
    P(f,2) = computePSNR(Iref, F{f,2});
    clean_f = imfilter(Iref, filt(f).h, 'replicate');
    gF = sobelMag(clean_f);
    fprintf('%-24s %10.2f %15.2f %18.2f %16.1f\n', filt(f).name, P(f,1), P(f,2), ...
            computePSNR(Iref, clean_f), 100*mean(gF(edgeMask))/mean(gRef(edgeMask)));
end
[~, b1] = max(P(:,1)); [~, b2] = max(P(:,2));
fprintf('Highest PSNR - Gaussian noise: %s (%.2f dB), salt & pepper: %s (%.2f dB)\n', ...
        filt(b1).name, P(b1,1), filt(b2).name, P(b2,2));

% check against MATLAB psnr() if it exists
try
    fprintf('Check: built-in psnr(noisy Gaussian) = %.2f dB\n', psnr(N.gauss, Iref));
catch
end
