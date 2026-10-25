%% EE3005 Image Processing - Assignment 2  -  run everything
%  Runs Task 1 ... Task 8. Figures go to the "figures" folder and the full
%  Command Window output is saved to "results/console_output.txt".
%  Every task file can also be run on its own (open it and press Run).
%
%  Needs: MATLAB + Image Processing Toolbox (no other toolbox).
thisDir = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
cd(thisDir);
addpath(fullfile(thisDir, 'helpers'));
P = getPaths();
logFile = fullfile(P.results, 'console_output.txt');
if exist(logFile, 'file'), delete(logFile); end
diary(logFile);
fprintf('EE3005 Assignment 2 - run started %s\n', datestr(now));

Task1_AHE_CLAHE;
Task2_Conv_Corr_Padding;
Task3_Separable;
Task4_LinearFiltering;
Task5_NonLinearFiltering;
Task6_UnsharpMasking;
Task7_EdgeDetection;
Task8_Pipeline;

fprintf('\nAll tasks finished. Figures: %s\n', getfield(getPaths(), 'figures'));
diary off;
