%% EE3005 Image Processing - Assignment 2
%  Task 5: Non-linear filtering and impulse noise removal
%  Same reference image and noisy images as Task 4 (makeNoisyImages)
clear; close all; clc;
thisDir = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'helpers'));
fprintf('\n================ TASK 5: NON-LINEAR FILTERING ================\n');

Iref = readGray('chest_xray.png');
N = makeNoisyImages(Iref);
roiZoom = {150:270, 120:240};
Q = 1.5;                                   % order used in the lecture example (Fig. 5.8)

% check of the geometric mean example from the lecture slide
blk = [5 16 22; 6 3 18; 12 3 15];
fprintf('Lecture example: (5*16*22*6*3*18*12*3*15)^(1/9) = %.2f\n', prod(blk(:))^(1/9));
gm = geometricMeanFilter(blk/255, 3) * 255;
fprintf('geometricMeanFilter on the same 3x3 block, centre value = %.2f\n\n', gm(2,2));

%% 5.1 Median filter vs arithmetic mean on salt & pepper noise
med3 = medfilt2(N.sp, [3 3], 'symmetric');   % 'symmetric' avoids black borders
med5 = medfilt2(N.sp, [5 5], 'symmetric');
mean3 = imfilter(N.sp, fspecial('average', 3), 'replicate');
mean5 = imfilter(N.sp, fspecial('average', 5), 'replicate');

fig = newFig('Task 5.1 - median vs mean', [40 40 1400 800]);
im51 = {N.sp, mean3, mean5, med3, med5, Iref};
nm51 = {'Salt & pepper noisy', 'Mean 3x3', 'Mean 5x5', 'Median 3x3', 'Median 5x5', 'Reference'};
for k = 1:6
    subplot(2,6,k);   imshow(im51{k}); title(nm51{k});
    subplot(2,6,k+6); imshow(im51{k}(roiZoom{:})); title('zoom');
end
saveFig(fig, 'Fig5_1_median_vs_mean');
fprintf('5.1 PSNR on salt & pepper: mean3 %.2f, mean5 %.2f, median3 %.2f, median5 %.2f dB\n', ...
        computePSNR(Iref, mean3), computePSNR(Iref, mean5), computePSNR(Iref, med3), computePSNR(Iref, med5));
%% 5.2 / 5.3 / 5.4  All filters on all four noise types
noisy = {N.gauss, N.sp, N.salt, N.pepper};
nName = {'Gaussian', 'Salt&pepper', 'Salt only', 'Pepper only'};
fName = {'Arithmetic mean 3x3', 'Median 3x3', 'Median 5x5', 'Geometric mean 3x3', ...
         'Harmonic mean 3x3', sprintf('Contra-harm. Q=+%.1f', Q), sprintf('Contra-harm. Q=-%.1f', Q)};
filters = {@(x) imfilter(x, fspecial('average', 3), 'replicate'), ...
           @(x) medfilt2(x, [3 3], 'symmetric'), ...
           @(x) medfilt2(x, [5 5], 'symmetric'), ...
           @(x) geometricMeanFilter(x, 3), ...
           @(x) harmonicMeanFilter(x, 3), ...
           @(x) contraharmonicFilter(x, 3,  Q), ...
           @(x) contraharmonicFilter(x, 3, -Q)};
R = cell(numel(filters), numel(noisy));
Ptab = zeros(numel(filters) + 1, numel(noisy));
for n = 1:numel(noisy)
    Ptab(1, n) = computePSNR(Iref, noisy{n});
    for f = 1:numel(filters)
        R{f, n} = min(max(filters{f}(noisy{n}), 0), 1);
        Ptab(f+1, n) = computePSNR(Iref, R{f, n});
    end
end

% Gaussian noise: geometric vs arithmetic mean, whole image and bright part only
bright = Iref > 0.4;                              % excludes dark lungs / black border
nz = mean(N.gauss(:) == 0);
psnrMask = @(a) 10*log10(1/mean((Iref(bright) - a(bright)).^2));
fprintf('\n5.2 Gaussian noise: %.2f %% of noisy pixels were clipped to 0\n', 100*nz);
fprintf('    PSNR whole image : arithmetic %.2f dB, geometric %.2f dB, harmonic %.2f dB\n', ...
        Ptab(2,1), Ptab(5,1), Ptab(6,1));
fprintf('    PSNR where ref > 0.4 : arithmetic %.2f dB, geometric %.2f dB, harmonic %.2f dB\n', ...
        psnrMask(R{1,1}), psnrMask(R{4,1}), psnrMask(R{5,1}));
fprintf('    mean error (filtered - ref) where ref > 0.4 : arithmetic %+.4f, geometric %+.4f, harmonic %+.4f\n', ...
        mean(R{1,1}(bright) - Iref(bright)), mean(R{4,1}(bright) - Iref(bright)), mean(R{5,1}(bright) - Iref(bright)));

%% 5.2 figures: each mean variant on an appropriate image and on a bad case
fig = newFig('Task 5.2 - mean filter variants', [40 40 1400 900]);
% row 1: geometric mean (good: Gaussian noise, bad: pepper)
subplot(3,4,1);  imshow(N.gauss);   title('Gaussian noisy');
subplot(3,4,2);  imshow(R{4,1});    title('Geometric 3x3 on Gaussian');
subplot(3,4,3);  imshow(N.pepper);  title('Pepper noisy');
subplot(3,4,4);  imshow(R{4,4});    title('Geometric 3x3 on pepper');
% row 2: harmonic mean (good: salt, bad: pepper)
subplot(3,4,5);  imshow(N.salt);    title('Salt noisy');
subplot(3,4,6);  imshow(R{5,3});    title('Harmonic 3x3 on salt');
subplot(3,4,7);  imshow(N.pepper);  title('Pepper noisy');
subplot(3,4,8);  imshow(R{5,4});    title('Harmonic 3x3 on pepper');
% row 3: contra-harmonic (good: Q>0 on pepper, bad: on salt & pepper)
subplot(3,4,9);  imshow(N.pepper);  title('Pepper noisy');
subplot(3,4,10); imshow(R{6,4});    title(sprintf('CH Q=+%.1f on pepper', Q));
subplot(3,4,11); imshow(N.sp);      title('Salt & pepper noisy');
subplot(3,4,12); imshow(R{6,2});    title(sprintf('CH Q=+%.1f on salt & pepper', Q));
saveFig(fig, 'Fig5_2_mean_variants');

fig = newFig('Task 5.2 - zoom', [40 40 1400 450]);
zz = {Iref, R{4,1}, R{1,1}, R{5,3}, R{4,4}, R{5,4}};
zn = {'Reference', 'Geometric (Gaussian)', 'Arithmetic (Gaussian)', 'Harmonic (salt)', ...
      'Geometric (pepper)', 'Harmonic (pepper)'};
for k = 1:6
    subplot(1,6,k); imshow(zz{k}(roiZoom{:})); title(zn{k});
end
saveFig(fig, 'Fig5_3_mean_variants_zoom');

%% 5.3 Contra-harmonic: sign of Q on salt and on pepper
fig = newFig('Task 5.3 - contra-harmonic', [40 40 1300 850]);
subplot(2,3,1); imshow(N.pepper); title('Pepper noise, p = 0.1');
subplot(2,3,2); imshow(R{6,4});   title(sprintf('Q = +%.1f  (correct sign)', Q));
subplot(2,3,3); imshow(R{7,4});   title(sprintf('Q = -%.1f  (wrong sign)', Q));
subplot(2,3,4); imshow(N.salt);   title('Salt noise, p = 0.1');
subplot(2,3,5); imshow(R{6,3});   title(sprintf('Q = +%.1f  (wrong sign)', Q));
subplot(2,3,6); imshow(R{7,3});   title(sprintf('Q = -%.1f  (correct sign)', Q));
saveFig(fig, 'Fig5_4_contraharmonic_sign');

% more Q values on pepper / salt to see the trend
fprintf('\n5.3 Contra-harmonic PSNR [dB] for several Q (3x3)\n');
fprintf('%8s %12s %12s\n', 'Q', 'Pepper', 'Salt');
for q = [-3 -1.5 -1 0 1 1.5 3]
    fprintf('%8.1f %12.2f %12.2f\n', q, ...
        computePSNR(Iref, min(max(contraharmonicFilter(N.pepper, 3, q), 0), 1)), ...
        computePSNR(Iref, min(max(contraharmonicFilter(N.salt,   3, q), 0), 1)));
end
fprintf('(Q = 0 is the arithmetic mean, Q = -1 is the harmonic mean)\n');

% edge effect of Q: mean brightness change on the clean image
fprintf('Effect on the CLEAN image: mean intensity change, Q=+%.1f: %+.4f, Q=-%.1f: %+.4f\n', ...
        Q, mean(mean(contraharmonicFilter(Iref,3,Q) - Iref)), Q, mean(mean(contraharmonicFilter(Iref,3,-Q) - Iref)));

%% 5.4 PSNR comparison table
fprintf('\n5.4 PSNR [dB] against the clean reference (all windows use replicate/symmetric borders)\n');
fprintf('%-24s', 'Filter'); fprintf('%13s', nName{:}); fprintf('\n');
fprintf('%-24s', 'Noisy (no filter)'); fprintf('%13.2f', Ptab(1,:)); fprintf('\n');
for f = 1:numel(filters)
    fprintf('%-24s', fName{f}); fprintf('%13.2f', Ptab(f+1,:)); fprintf('\n');
end
[~, best] = max(Ptab(2:end, :), [], 1);
for n = 1:numel(noisy)
    fprintf('Best for %-12s: %s (%.2f dB)\n', nName{n}, fName{best(n)}, Ptab(best(n)+1, n));
end
