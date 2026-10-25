%% EE3005 Image Processing - Assignment 2
%  Task 6: Image sharpening with unsharp masking (and sharp masking)
%  Image: chest X-ray (ribs, vessels, heart border)
clear; close all; clc;
thisDir = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'helpers'));
fprintf('\n================ TASK 6: UNSHARP MASKING ================\n');

f = readGray('chest_xray.png');
roiZoom = {150:270, 120:240};          % right lung, vessels and ribs
roiFlat = {420:470, 100:170};          % flat liver region (noise / grain)
profRow = 250; profCols = 300:400;     % line crossing the left heart border

%% 6.1 Blurred image and mask
sigmaBlur = 2;
hBlur = gaussKernel(sigmaBlur);        % 13x13 Gaussian, sigma = 2
fBlur = imfilter(f, hBlur, 'replicate');
m = f - fBlur;                         % unsharp mask m = f - f_blur

%% 6.2 Several gains
gammas = [0.5 1 2 4];
g = cell(size(gammas));
for k = 1:numel(gammas)
    g{k} = f + gammas(k) * m;          % g = f + gamma*m (not clipped yet)
end

fig = newFig('Task 6.1 - blur and mask', [40 40 1300 450]);
subplot(1,3,1); imshow(f);     title('Original f');
subplot(1,3,2); imshow(fBlur); title(sprintf('Blurred, Gaussian sigma = %d', sigmaBlur));
subplot(1,3,3); imshow(m, []); title('Mask m = f - f_{blur} (scaled)');
saveFig(fig, 'Fig6_1_blur_mask');

fig = newFig('Task 6.2 - sharpened', [40 40 1400 800]);
for k = 1:numel(gammas)
    subplot(2,4,k);   imshow(min(max(g{k},0),1)); title(sprintf('\\gamma = %.1f', gammas(k)));
    subplot(2,4,k+4); imshow(min(max(g{k}(roiZoom{:}),0),1)); title(sprintf('zoom, \\gamma = %.1f', gammas(k)));
end
saveFig(fig, 'Fig6_2_sharpened_gammas');

% intensity profile across the heart border
fig = newFig('Task 6.2 - profile', [40 40 1000 450]);
plot(profCols, f(profRow, profCols), 'k', 'LineWidth', 2); hold on;
cols = [0.2 0.4 0.8; 0.1 0.6 0.2; 0.9 0.5 0.1; 0.8 0.1 0.1];
for k = 1:numel(gammas)
    plot(profCols, g{k}(profRow, profCols), 'Color', cols(k,:), 'LineWidth', 1.2);
end
hold off; grid on; xlim([profCols(1) profCols(end)]);
xlabel('column'); ylabel('intensity (before clipping)');
legend([{'original'}, arrayfun(@(x) sprintf('\\gamma = %.1f', x), gammas, 'UniformOutput', false)], ...
       'Location', 'best');
title(sprintf('Profile along row %d (heart border)', profRow));
saveFig(fig, 'Fig6_3_edge_profile');

% numbers for the discussion
gradF = sobelMag(f);
fprintf('6.2 Effect of gamma (blur sigma = %d)\n', sigmaBlur);
fprintf('%-18s %12s %14s %14s %12s\n', 'Image', 'Mean |grad|', 'Flat-ROI noise', ...
        'Out of [0,1] %', 'Overshoot');
fprintf('%-18s %12.4f %14.4f %14s %12s\n', 'Original', mean(gradF(:)), ...
        highpassStd(f, roiFlat{1}, roiFlat{2}), '0.00', '-');
lp = f(profRow, profCols); rangeF = max(lp) - min(lp);
for k = 1:numel(gammas)
    gk = g{k}; gc = min(max(gk, 0), 1); gr = sobelMag(gc);
    pk = gk(profRow, profCols);
    over = (max(pk) - max(lp) + min(lp) - min(pk)) / rangeF;   % halo size on the profile
    fprintf('%-18s %12.4f %14.4f %14.2f %11.0f%%\n', sprintf('gamma = %.1f', gammas(k)), ...
        mean(gr(:)), highpassStd(gc, roiFlat{1}, roiFlat{2}), ...
        100*mean(gk(:) < 0 | gk(:) > 1), 100*over);
end
fprintf('(Overshoot = extra height of the profile above/below the original edge, %% of edge height)\n');

%% 6.3 Unsharp masking as one linear kernel, and sharp masking
gSel = 1;                                           % selected gain
delta = zeros(size(hBlur)); delta((end+1)/2, (end+1)/2) = 1;
hUnsharp = (1 + gSel) * delta - gSel * hBlur;       % g = (1+gamma) f - gamma f_blur
gOne = imfilter(f, hUnsharp, 'replicate');
gTwo = f + gSel * m;
fprintf('\n6.3 One-kernel unsharp masking vs two-step: max difference = %.2e\n', max(abs(gOne(:) - gTwo(:))));

% sharp masking: high-pass detail taken directly (no blur step)
% high-pass = -Laplacian, kernel [0 -1 0; -1 4 -1; 0 -1 0]; g = f + a * highpass
hHP = [0 -1 0; -1 4 -1; 0 -1 0];
hp = imfilter(f, hHP, 'replicate');
aList = [0.5 1];
gs = cell(size(aList));
for k = 1:numel(aList)
    gs{k} = f + aList(k) * hp;                     % a = 1 gives kernel [0 -1 0; -1 5 -1; 0 -1 0]
end
fig = newFig('Task 6.3 - unsharp vs sharp masking', [40 40 1400 800]);
cmp  = {f, min(max(g{2},0),1), min(max(gs{1},0),1), min(max(gs{2},0),1)};
cmpN = {'Original', 'Unsharp, \gamma = 1', 'Sharp masking, a = 0.5', 'Sharp masking, a = 1'};
for k = 1:4
    subplot(2,4,k);   imshow(cmp{k}); title(cmpN{k});
    subplot(2,4,k+4); imshow(cmp{k}(roiZoom{:})); title('zoom');
end
saveFig(fig, 'Fig6_4_unsharp_vs_sharp');

fprintf('\nUnsharp vs sharp masking\n');
fprintf('%-26s %12s %14s %14s\n', 'Image', 'Mean |grad|', 'Flat-ROI noise', 'Out of [0,1] %');
items = {g{2}, gs{1}, gs{2}}; itemsN = {'Unsharp gamma=1 (s=2)', 'Sharp masking a=0.5', 'Sharp masking a=1'};
for k = 1:3
    gc = min(max(items{k}, 0), 1); gr = sobelMag(gc);
    fprintf('%-26s %12.4f %14.4f %14.2f\n', itemsN{k}, mean(gr(:)), ...
            highpassStd(gc, roiFlat{1}, roiFlat{2}), 100*mean(items{k}(:) < 0 | items{k}(:) > 1));
end
% how much low-frequency (large-scale) content each detail layer contains
lowM  = imfilter(m,  gaussKernel(4), 'replicate');
lowHP = imfilter(hp, gaussKernel(4), 'replicate');
fprintf('Share of detail energy at large scale (after sigma=4 blur): unsharp mask %.1f %%, high-pass %.1f %%\n', ...
        100*sum(lowM(:).^2)/sum(m(:).^2), 100*sum(lowHP(:).^2)/sum(hp(:).^2));
