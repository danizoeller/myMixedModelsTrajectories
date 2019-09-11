function [fh,residuals] = plotResiduals(model,figurehandle)
% plotResiduals - residuals norm plot
%
% Syntax:  [fh,residuals] = plotResiduals(model,figurehandle)
%
% Inputs:
%    model - structure containing estimated mixed-effect model, see output
%       of estimateModel() for more details
%    plotOpts - struct with plotting color and the axis label and legend 
%       information, if empty, the variable names will be used as labels
%       [.plotType] = 'full' - always plot full model, even if interation or 
%                       intercept are not significantly different
%                'redInter' - plot reduced model without interaction if 
%                       interaction is not significant
%                'redGrp' - plot reduced model without group effect if 
%                       group effect is not significant 
%       [.plotCol] - cellstring with plotting colors
%    figurehandle - handle to figure where to plot the graph, if emty or
%       not existent a new figure is created
%
% Outputs:
%    fh - handle of plotted figure
%    plottedModel - string indicating which model was plotted:
%           'full','noGroup' or 'noInteraction'
%
% Example: 
%

% Author: Daniela Zoeller
% Medical Image Processing Lab, EPFL/UniGe
% Developmental Imaging and Psychopathology Lab, UniGe
% v1.0 3.8.2016 DZ - initial version based on initial code of Kadir Mutlu

%% set up defaults
if nargin < 2 || isempty(figurehandle)
    fh=figure;
else
    fh = figurehandle;
end

nCov = length(model.beta)-2*(model.order+1);

%% residuals
pred=model.designMatrix*model.beta;
residuals = pred - model.input.data;

%% normplot
figure(fh)
normplot(residuals)





