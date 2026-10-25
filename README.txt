EE3005 Image Processing - Assignment 2 : MATLAB code
=====================================================

HOW TO RUN
1. Unzip the folder anywhere (keep the folder structure).
2. Open MATLAB, go to this folder, open START_HERE.m and press Run.
   -> runs Task 1..8, saves all figures to  figures\
   -> saves the whole Command Window output to  results\console_output.txt
3. Or open any TaskX_....m file and press Run to run only that task.
   (each task file adds the helpers folder to the path by itself)

Requirement: MATLAB with Image Processing Toolbox. No other toolbox.

FILES
Task1_AHE_CLAHE.m           global HE, AHE (adapthisteq ClipLimit=1), CLAHE parameters
Task2_Conv_Corr_Padding.m   hand calculation check (Figure 1a), 4 padding types
Task3_Separable.m           SVD separability, 2-D vs 1-D+1-D, operation count, timing
Task4_LinearFiltering.m     Gaussian / salt&pepper noise, mean + Gaussian filters, PSNR
Task5_NonLinearFiltering.m  median, geometric, harmonic, contra-harmonic (Q=+-1.5)
Task6_UnsharpMasking.m      unsharp masking (gamma 0.5,1,2,4), sharp masking
Task7_EdgeDetection.m       Sobel, Laplacian, LoG, DoG, noise sensitivity
Task8_Pipeline.m            median -> CLAHE -> unsharp, plus two other orders
helpers\                    small functions used by the tasks (PSNR, filters, ...)
images\chest_xray.png       NIH ChestX-ray14 image 00000001_000 (512x512, 8-bit)
images\brain_mri.png        axial slice z=80 of the MNI ICBM152 2009a T1 template

NOTE
Noise is generated with fixed random seeds (rng), so every run gives the
same noisy images and the same numbers.
