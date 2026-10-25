%% EE3005 Image Processing - Assignment 2
%  Task 8: Integrated pipeline  (noise reduction -> CLAHE -> unsharp masking)
%  Test image: brain MRI made low-contrast and noisy on purpose
clear; close all; clc;
thisDir = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'helpers'));
fprintf('\n================ TASK 8: PIPELINE ================\n');

Iclean = readGray('brain_mri.png');

%% Degraded input (artificial)
% 1) contrast compressed to [0.30, 0.60]
% 2) Gaussian noise, sigma = 0.01 (variance 1e-4)
% 3) salt & pepper noise, density 0.02
Ilow = 0.30 + 0.30 * Iclean;
rng(10); Ideg = Ilow + 0.01 * randn(size(Ilow));
Ideg = min(max(Ideg, 0), 1);
rng(11); Ideg = imnoise(Ideg, 'salt & pepper', 0.02);

%% Pipeline A (proposed): median -> CLAHE -> unsharp masking
S1 = medfilt2(Ideg, [3 3], 'symmetric');                        % noise reduction (Task 5)
S2 = adapthisteq(S1, 'NumTiles', [8 8], 'ClipLimit', 0.01);     % contrast (Task 1)
gammaU = 0.7; sigmaU = 1.5;
S3 = S2 + gammaU * (S2 - imfilter(S2, gaussKernel(sigmaU), 'replicate'));   % sharpening (Task 6)
S3 = min(max(S3, 0), 1);

%% Alternative order B: CLAHE -> median -> unsharp masking
B1 = adapthisteq(Ideg, 'NumTiles', [8 8], 'ClipLimit', 0.01);
B2 = medfilt2(B1, [3 3], 'symmetric');
B3 = min(max(B2 + gammaU * (B2 - imfilter(B2, gaussKernel(sigmaU), 'replicate')), 0), 1);

%% Alternative order C: median -> unsharp masking -> CLAHE
C1 = S1;
C2 = min(max(C1 + gammaU * (C1 - imfilter(C1, gaussKernel(sigmaU), 'replicate')), 0), 1);
C3 = adapthisteq(C2, 'NumTiles', [8 8], 'ClipLimit', 0.01);

%% Figures
fig = newFig('Task 8 - pipeline stages', [40 40 1400 800]);
st = {Iclean, Ideg, S1, S2, S3};
sn = {'Clean reference (not used)', 'Input: low contrast + noise', '1) Median 3x3', ...
      '2) + CLAHE (0.01, 8x8)', sprintf('3) + Unsharp (\\sigma=%.1f, \\gamma=%.1f)', sigmaU, gammaU)};
for k = 1:5
    subplot(2,5,k);   imshow(st{k}, [0 1]); title(sn{k});
    subplot(2,5,k+5); showHist(st{k}, 'histogram');
end
saveFig(fig, 'Fig8_1_pipeline_stages');

fig = newFig('Task 8 - order comparison', [40 40 1400 450]);
od = {Ideg, S3, B3, C3};
on = {'Input', 'A: median-CLAHE-unsharp', 'B: CLAHE-median-unsharp', 'C: median-unsharp-CLAHE'};
for k = 1:4
    subplot(1,4,k); imshow(od{k}, [0 1]); title(on{k});
end
saveFig(fig, 'Fig8_2_order_comparison');

roi = {60:140, 50:150};
fig = newFig('Task 8 - zoom', [40 40 1400 400]);
for k = 1:4
    subplot(1,4,k); imshow(od{k}(roi{:}), [0 1]); title([on{k} ' (zoom)']);
end
saveFig(fig, 'Fig8_3_zoom');

%% Evaluation
% tissue masks from the CLEAN image (only used for measuring)
se = ones(3);
WM  = imerode(Iclean >= 0.88, se);                     % white matter
GM  = imerode(Iclean >= 0.62 & Iclean <= 0.78, se);    % grey matter
CSF = imerode(Iclean > 0.10 & Iclean < 0.40, se);      % ventricles / CSF
BG  = imerode(Iclean < 0.02, ones(11));                % background
gC = sobelMag(Iclean); sv = sort(gC(:)); edgeM = gC >= sv(round(0.95*numel(sv)));

fprintf('Masks: WM %d px, GM %d px, CSF %d px, background %d px\n', ...
        nnz(WM), nnz(GM), nnz(CSF), nnz(BG));
fprintf('%-30s %9s %9s %9s %9s %10s %9s\n', 'Image', 'WM-GM', 'CNR', 'WM-CSF', 'CNR', 'Noise WM', 'Noise BG');
allI = {Iclean, Ideg, S1, S2, S3, B3, C3};
allN = {'Clean reference', 'Degraded input', 'A1 median', 'A2 median+CLAHE', ...
        'A3 final (proposed)', 'B CLAHE-median-unsharp', 'C median-unsharp-CLAHE'};
for k = 1:numel(allI)
    X = allI{k};
    nWM = std(X(WM));
    d1 = mean(X(WM)) - mean(X(GM)); d2 = mean(X(WM)) - mean(X(CSF));
    fprintf('%-30s %9.3f %9.1f %9.3f %9.1f %10.4f %9.4f\n', allN{k}, d1, d1/nWM, d2, d2/nWM, nWM, std(X(BG)));
end
fprintf('(CNR = mean difference / std inside white matter; noise BG = std in the black background)\n');

fprintf('\nEdge strength on the 5 %% strongest true edges (mean Sobel magnitude):\n');
for k = 1:numel(allI)
    gk = sobelMag(allI{k});
    fprintf('   %-30s %.3f\n', allN{k}, mean(gk(edgeM)));
end
fprintf('Impulses left (pixels exactly 0 or 1 inside the head): input %d, A %d, B %d, C %d\n', ...
        nnz((Ideg == 0 | Ideg == 1) & ~BG), nnz((S1 == 0 | S1 == 1) & ~BG), ...
        nnz((B2 == 0 | B2 == 1) & ~BG), nnz((C1 == 0 | C1 == 1) & ~BG));
