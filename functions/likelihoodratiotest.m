function [h,p,dof,lratio] = likelihoodratiotest(llogModel,llogNull,dof,alpha)
% likelihoodratiotest - function for likelihood ratio test calculation
%
% Syntax:  [h,p] = likelihoodratiotest(llogNull,llogModel,dof,[alpha])
%
% Inputs:
%    llogModel - log-likelihood of full model
%    llogNull - log-likelihood of reduced model
%    dof - degrees of freedom of Chi-Square distribution: difference of
%       parameter numbers between the two models
%    [alpha] - alpha for significance testing, default = 0.05
%
% Outputs:
%    h - 1: null hypothesis rejected, i.e. full model is significantly
%       better than the reduced model
%    p - p-value of test statistics
%
% Example: 
%

% Author: Daniela Zoeller
% Medical Image Processing Lab, EPFL/UniGe
% Developmental Imaging and Psychopathology Lab, UniGe
% v1.0 3.8.2016 DZ - initial version based on initial code of Kadir Mutlu

if nargin < 4
    alpha=0.05;
end

lratio = 2*(llogModel-llogNull);
p = 1-chi2cdf(lratio,dof);
h = (p <= alpha);