%% EE3005 Image Processing - Assignment 2
%  Task 1: Global HE, Adaptive Histogram Equalization (AHE) and CLAHE
%  Image  : chest X-ray (images/chest_xray.png, 512x512, 8-bit)
clear; close all; clc;
thisDir = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'helpers'));
fprintf('\n================ TASK 1: HE / AHE / CLAHE ================\n');

I8 = im2uint8(readGray('chest_xray.png'));

% Regions of interest (row range, column range) used in the analysis
roiLung  = {120:300,  90:210};   % right lung (left side of the image)
roiHeart = {200:350, 220:330};   % heart / mediastinum (dense, bright)
roiFlat  = {420:470, 100:170};   % liver area below diaphragm (nearly flat)

%% 1.1 Original image and histogram
fig = newFig('Task 1.1 - original', [40 40 1150 450]);
subplot(1,2,1); imshow(I8); title('Original chest X-ray'); hold on;
rectangle('Position', [roiLung{2}(1)  roiLung{1}(1)  numel(roiLung{2})  numel(roiLung{1})],  'EdgeColor', 'y', 'LineWidth', 1.2);
rectangle('Position', [roiHeart{2}(1) roiHeart{1}(1) numel(roiHeart{2}) numel(roiHeart{1})], 'EdgeColor', 'c', 'LineWidth', 1.2);
rectangle('Position', [roiFlat{2}(1)  roiFlat{1}(1)  numel(roiFlat{2})  numel(roiFlat{1})],  'EdgeColor', 'm', 'LineWidth', 1.2);
hold off;
subplot(1,2,2); showHist(I8, 'Histogram of the original image');
saveFig(fig, 'Fig1_1_original_histogram');

sv = sort(double(I8(:)));                 % percentiles without extra toolboxes
pct = sv(max(1, round([1 5 50 95 99]/100*numel(sv))))';
fprintf('1.1 Original: min=%d max=%d, percentiles 1/5/50/95/99 = %s\n', ...
        min(I8(:)), max(I8(:)), mat2str(pct));
fprintf('    mean grey level: lung ROI = %.1f, heart ROI = %.1f, flat ROI = %.1f\n', ...
        mean(mean(double(I8(roiLung{:})))), mean(mean(double(I8(roiHeart{:})))), ...
        mean(mean(double(I8(roiFlat{:})))));
fprintf('    std inside ROI: lung = %.1f, heart = %.1f (grey levels)\n', ...
        std(double(reshape(I8(roiLung{:}),[],1))), std(double(reshape(I8(roiHeart{:}),[],1))));

%% 1.2 Global HE and AHE
J_he  = im2uint8(histeq(I8));                                       % MATLAB default: 64 output levels
J_ahe = adapthisteq(I8, 'NumTiles', [8 8], 'ClipLimit', 1); % ClipLimit = 1 -> no clipping = plain AHE

fig = newFig('Task 1.2 - HE and AHE', [40 40 1300 800]);
subplot(2,3,1); imshow(I8);    title('Original');
subplot(2,3,2); imshow(J_he);  title('Global HE (histeq)');
subplot(2,3,3); imshow(J_ahe); title('AHE (8x8 tiles, no clip limit)');
subplot(2,3,4); showHist(I8,    'Original');
subplot(2,3,5); showHist(J_he,  'Global HE');
subplot(2,3,6); showHist(J_ahe, 'AHE');
saveFig(fig, 'Fig1_2_HE_AHE');

% zoomed regions: lung and heart
fig = newFig('Task 1.2 - zoom', [40 40 1150 750]);
imgs = {I8, J_he, J_ahe}; names = {'Original', 'Global HE', 'AHE'};
for k = 1:3
    subplot(2,3,k);   imshow(imgs{k}(roiLung{:}));  title([names{k} ' - lung ROI']);
    subplot(2,3,k+3); imshow(imgs{k}(roiHeart{:})); title([names{k} ' - heart ROI']);
end
saveFig(fig, 'Fig1_3_zoom_HE_AHE');

%% 1.3 CLAHE with different parameters
% effect of the clip limit (tiles fixed to 8x8)
clipList = [0.005 0.01 0.03];
% effect of the tile size (clip limit fixed to 0.01)
tileList = {[4 4], [16 16]};

J_clahe = {}; labels = {};
for c = clipList
    J_clahe{end+1} = adapthisteq(I8, 'NumTiles', [8 8], 'ClipLimit', c); %#ok<SAGROW>
    labels{end+1}  = sprintf('CLAHE clip=%.3f, 8x8', c);                %#ok<SAGROW>
end
for t = 1:numel(tileList)
    J_clahe{end+1} = adapthisteq(I8, 'NumTiles', tileList{t}, 'ClipLimit', 0.01); %#ok<SAGROW>
    labels{end+1}  = sprintf('CLAHE clip=0.010, %dx%d', tileList{t});             %#ok<SAGROW>
end

fig = newFig('Task 1.3 - CLAHE clip limit', [40 40 1300 800]);
for k = 1:3
    subplot(2,3,k);   imshow(J_clahe{k}); title(labels{k});
    subplot(2,3,k+3); showHist(J_clahe{k}, labels{k});
end
saveFig(fig, 'Fig1_4_CLAHE_cliplimit');

fig = newFig('Task 1.3 - CLAHE tile size', [40 40 1300 800]);
order = [4 2 5];                        % 4x4, 8x8, 16x16 (all clip 0.01)
for k = 1:3
    subplot(2,3,k);   imshow(J_clahe{order(k)}); title(labels{order(k)});
    subplot(2,3,k+3); imshow(J_clahe{order(k)}(roiLung{:})); title('lung ROI');
end
saveFig(fig, 'Fig1_5_CLAHE_tilesize');

%% 1.4 Comparison table
allImgs   = [{I8, J_he, J_ahe}, J_clahe];
allLabels = [{'Original', 'Global HE', 'AHE (no clip)'}, labels];
fprintf('\n1.4 Comparison (grey levels 0-255)\n');
fprintf('%-26s %8s %10s %9s %10s %10s %10s\n', 'Method', 'Mean', 'Global std', ...
        'Entropy', 'LocStd lung', 'LocStd heart', 'Flat noise');
for k = 1:numel(allImgs)
    D = double(allImgs{k});
    locStd = stdfilt(D, ones(9));          % local contrast in 9x9 windows
    fprintf('%-26s %8.1f %10.1f %9.2f %10.2f %10.2f %10.2f\n', allLabels{k}, mean(D(:)), ...
        std(D(:)), entropy(allImgs{k}), mean(mean(locStd(roiLung{:}))), ...
        mean(mean(locStd(roiHeart{:}))), highpassStd(D, roiFlat{1}, roiFlat{2}));
end
fprintf('(LocStd = mean 9x9 local standard deviation inside the ROI,\n');
fprintf(' Flat noise = std of fine-scale fluctuation in the flat liver ROI)\n');

% final side-by-side comparison
best = 2;                                 % CLAHE clip 0.01, 8x8 (chosen setting)
fig = newFig('Task 1.4 - comparison', [40 40 1400 700]);
cmp = {I8, J_he, J_ahe, J_clahe{best}}; cmpN = {'Original', 'Global HE', 'AHE', labels{best}};
for k = 1:4
    subplot(2,4,k);   imshow(cmp{k}); title(cmpN{k});
    subplot(2,4,k+4); imshow(cmp{k}(roiFlat{1}(1)-60:roiFlat{1}(end), roiFlat{2}(1)-20:roiFlat{2}(end)+60));
    title('lower-left region (liver / diaphragm)');
end
saveFig(fig, 'Fig1_6_comparison');
