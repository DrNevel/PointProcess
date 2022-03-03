function [SD1, SD2, SD_ratio, SD_prod, Var, coeff, x, y]=get_poincare(hrv)
% Computes SD1 and SD2 [ms] as well as providing an ellipse E_c centered in the x-y
% mean of the poincare plot. In addition it returns the variances resulted
% from the application of pca in the column vector Var [s2], in which the first
% row refers to SD2 and the second to SD1. coeff are the eigenvalues of the
% pca. SD1 and SD2 are obtained through the application of Principal
% Component Analysis of the scatter plot.

% SD1 - SD2
x=hrv(1:end-1);
y=hrv(2:end);
mn=[mean(x) mean(y)];
X=[x-mean(x); y-mean(y)]';
[coeff, ~, latent]=pca(X);
std_c=sqrt(latent);
SD2=std_c(1)*1000; % Horizontal Direction
SD1=std_c(2)*1000; % Vertical Direction
Var=latent;

% SD1/SD2
SD_ratio = SD1/SD2;

% SD1*SD2
SD_prod = SD1*SD2;


end