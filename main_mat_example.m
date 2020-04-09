%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example main script to use the mixed effect models toolbox on data stored
% in an examplary mat file
%
% This script gives an example for using the code to fit mixed models
% trajectories with input data from an mat file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close and clear everything
clc
clear all
close all

% set up path to all necessary functions
addpath(genpath('./functions'));
startup; % set up some options for nicer plots

%% ------------------------------------------------------------------------
% set up all necessary otions here (input file, output directory, etc.)
% ------------------------------------------------------------------------

% -----------------------------
% input data options
% -----------------------------
inputDataFile='exampleData.mat'; % filename and path of excel file with input data
load(inputDataFile);

input.subjID=within_subj; % subject IDs
input.age=age; % age
input.grouping=diagnosis;% column (with 0 and 1 values) you want to use for grouping: $
                          % input.grouping=[]; if you have only 1 group
                          % 1 column if you have 2 groups
                          % 2 colums if you have 3 groups (in the first column  you have 1 for everyone in the group 1 and 0 everywhere else,
                          % and in the second column you have 1 for everyone is in the group 2 and 0 everywhere else)             
input.data=lh_thickness; % data to fit (thickness, volume, behavior, ...)
input.cov=sex; % model covariates (here, only sex is included as covariate)

% -----------------------------
% model estimation options
% -----------------------------
opts.orders=0:3; % model orders to check for (0=constant, 1=linear, 2=quadratic, 3=cubic, etc. ...)
opts.mType = 'slope'; % 'intercept' for random intercept, 'slope' for random slope (recommended)

% vertex IDs (colums of lh_thickness to analyze)
opts.vertID=[2000,3600,4000];

% -----------------------------
% model plotting options
% -----------------------------
outDir = fullfile('./results_mat_fdrcorr'); % directory where to store the result table and plots
saveResults = 2; % do you want to save the results? 
                 % 0=No; 1=Yes, but only the table; 2=Yes, both the table and the plots
plotOpts.legTxt = {'HC','Pat'}; % legend: names of your groups (here 0='HC' and 1='Pat')
plotOpts.xLabel = 'age'; % label for x-axis
plotOpts.yLabel = 'cortical thickness'; % label for y-axis
plotOpts.plotCI = 1; % do you want to plot confidence intervals? 0=No; 1=Yes
plotOpts.plotType = 'redInter'; % which models do you want to plot?
        % 'full' - always plot full model, even if interation or intercept are not significantly different
        % 'redInter' - plot reduced model without interaction if interaction is not significant
        % 'redGrp' - plot reduced model without group effect if group effect is not significant 

% you can also try some different colors by uncommenting the next two lines :-)
% colors = cbrewer('qual', 'Set2', 3);
% plotOpts.plotCol = {colors(1,:),colors(2,:),colors(3,:)}; 


%% ------------------------------------------------------------------------
% execute the model estimation and plot/save results
% ------------------------------------------------------------------------

% -----------------------------
% run model fitting
% -----------------------------
% a few more options
opts.modelNames = textscan(num2str(opts.vertID),'%s'); 
opts.modelNames=opts.modelNames{1};
opts.alpha=0.05; % significance level for group and interaction effects
opts.figPosition=[440   488   525   310]; % change the size of your figure (take a figure, change the size to whatever you want, type get(gcf,'Position') and copy the resulting output here


% fit models
outModelVect = fitOptModel(input,opts);

% correct for multiple comparisons using FDR
outModelVect_corr = fdr_correct(outModelVect,opts.alpha);

% -----------------------------
% plot and save models
% -----------------------------
plotOpts.nCov=size(input.cov,2);
plotModelsAndSaveResults(outModelVect_corr,plotOpts,saveResults,outDir);

% -----------------------------
% calculation of the size effect
% -----------------------------
% Reporting the effect size for the groups and for the interaction between groups and age

effectSizeGroup=GroupCalculationEffect(outModelVect);
%table reporting the group (if variable group is included) size effect for each model (each response variable)
effectSizeInter=InterCalculationEffect(outModelVect);
%table reporting the interaction (if variable interaction is included) size effect for each model (each response variable)


