%% EE3005 Image Processing - Assignment 2
%  Task 3: Separable filtering and computational efficiency
clear; close all; clc;
thisDir = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'helpers'));
fprintf('\n================ TASK 3: SEPARABLE FILTERING ================\n');

%% 3.1 Separability check of the kernel in Figure 1(b)
K = [1 2 1; 2 4 2; 1 2 1] / 16;
h = [1; 2; 1] / 4;        % column vector  -> filters along the vertical direction
v = [1 2 1] / 4;          % row vector     -> filters along the horizontal direction

s = svd(K);               % lecture: separable if only sigma_0 is non-zero
fprintf('Singular values of K: %s\n', mat2str(s', 4));
fprintf('rank(K) = %d\n', rank(K));
Kouter = h * v;           % outer product (3x1)*(1x3) = 3x3
fprintf('h*v * 16 =\n'); disp(Kouter * 16);
fprintf('max|h*v - K| = %g  -> identical: %d\n', max(abs(Kouter(:) - K(:))), isequal(Kouter, K));

% the SVD also gives the 1-D factors directly: K = s1 * u1 * v1'
[U, S, V] = svd(K);
hS = U(:,1) * sqrt(S(1,1));  vS = V(:,1)' * sqrt(S(1,1));
fprintf('1-D factors from SVD (sign may be negative): h = %s, v = %s\n', ...
        mat2str(hS', 4), mat2str(vS, 4));

%% 3.2 2-D filtering vs two 1-D filters
I = readGray('chest_xray.png');               % double in [0,1]
out2D  = imfilter(I, K, 'replicate');          % one 3x3 pass
tmp    = imfilter(I, v, 'replicate');          % horizontal pass (1x3)
outSep = imfilter(tmp, h, 'replicate');        % vertical pass   (3x1)
D = abs(out2D - outSep);
fprintf('\n3.2 double image: max |2-D - separable| = %.3e, mean = %.3e\n', max(D(:)), mean(D(:)));

% same test with uint8 data (the result is rounded after each pass)
I8 = im2uint8(I);
o2D8  = imfilter(I8, K, 'replicate');
oSep8 = imfilter(imfilter(I8, v, 'replicate'), h, 'replicate');
D8 = abs(double(o2D8) - double(oSep8));
fprintf('    uint8 image : max |2-D - separable| = %d grey level(s), %.2f %% of pixels differ\n', ...
        max(D8(:)), 100*mean(D8(:) > 0));

fig = newFig('Task 3.2 - 2D vs separable', [40 40 1400 420]);
subplot(1,4,1); imshow(I);      title('Original');
subplot(1,4,2); imshow(out2D);  title('2-D kernel K (3x3)');
subplot(1,4,3); imshow(outSep); title('h then v (1-D, 1-D)');
subplot(1,4,4); imshow(D, []);  title(sprintf('|difference|, max = %.1e', max(D(:)))); colorbar;
saveFig(fig, 'Fig3_1_separable');

%% 3.3 Number of multiplications per pixel and timing
fprintf('\n3.3 Multiplications per output pixel\n');
fprintf('%6s %12s %16s %10s\n', 'K', '2-D (K^2)', 'separable (2K)', 'ratio');
for Ks = [3 5 7 9 15 31]
    fprintf('%6d %12d %16d %10.1f\n', Ks, Ks^2, 2*Ks, Ks^2/(2*Ks));
end

% timing with conv2 (conv2 does exactly what we ask, it does not
% decide on its own to use separable filtering)
nRep = 50;
g31 = fspecial('gaussian', [31 1], 5);                % 1-D Gaussian, length 31
kernels = {h, g31};                                    % 1-D column factors
for k = 1:2
    g1 = kernels{k};
    G2 = g1 * g1';                                     % equivalent 2-D kernel
    a = conv2(I, G2, 'same'); b = conv2(g1, g1', I, 'same');   % warm-up
    tic; for r = 1:nRep, a = conv2(I, G2, 'same');         end; t2 = toc/nRep;
    tic; for r = 1:nRep, b = conv2(g1, g1', I, 'same');    end; ts = toc/nRep;
    fprintf('Kernel %2dx%-2d : 2-D %.2f ms, separable %.2f ms, speed-up %.1fx, max diff %.1e\n', ...
            numel(g1), numel(g1), 1000*t2, 1000*ts, t2/ts, max(abs(a(:) - b(:))));
end
fprintf('(timing depends on the computer; MATLAB conv2(u,v,A) runs the two 1-D passes)\n');
