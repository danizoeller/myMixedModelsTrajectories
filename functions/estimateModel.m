function outModel = estimateModel(subjID, grouping, age, cov, data, opts)
% estimateModel - function to fit a random intersect model of given order
%
% Syntax:  model = estimateModel(subjID, gouping, age, data, mOrder)
%
% Inputs:
%    subjID - vector of subject IDs (double, #sub x 1)
%    age - vector of ages (double, #sub x 1)
%    grouping - 1 vector (if 2 groups) or 2 vectors (if 3 groups) with grouping information 
%    cov - matrix of covariates (#sub < #cov), empty if no covariates
%    data - vector with data to model (double, #sub x 1)
%    opts - structure with options for model estimation
%           .mOrder - order of model to estimate
%           [.groupEffect] - true if grouping information should be 
%               included in the model, default: true
%           [.interEffect] - true if age*grouping interaction should be
%               considered, default: true
%           [.mType] - 'intersect' (default), random intersect model
%                      'slope' for random slope model
%
% Outputs:
%    outModel - structure with fitted model
%           .input - estimation inputs
%           .order - order of resulting model (=mOrder)
%           .designMatrix - design matrix of estimated model
%           .designVars - names of columns of the design matrix
%           .beta - vector of estimated fixed effects (betas)
%           .randCov - estimated covariance matrix for random effects
%           .stats  - structure with model fitting statistics (see nlmefit
%                   for details)
%
% Example: 
%

% Author: Daniela Zoeller
% Medical Image Processing Lab, EPFL/UniGe
% Developmental Imaging and Psychopathology Lab, UniGe
% v1.0 3.8.2016 DZ - initial version based on initial code of Kadir Mutlu
%%

%% set up defaults
if ~isfield(opts,'groupEffect')
    opts.groupEffect=1;
end
if ~isfield(opts,'interEffect')
    opts.interEffect=1;
end
if ~isfield(opts,'mType')
    opts.mType='intercept';
end
if ~opts.groupEffect
    opts.interEffect=0;
end

%% get constants
nObs = length(subjID);
mOrder = opts.mOrder;
% calculate the number of groups
if size(grouping,2)==0 
    groups=1;
else
groups=size(grouping,2)+1;
end

if groups==2 && ~isequal(unique(grouping),[0 1]')
    error('error occurred in the specification of the groups')
end
if groups==3 && (~isequal(unique(grouping(:,1)),[0 1]') || ~isequal(unique(grouping(:,2)),[0 1]'))
    error('error occurred in the specification of the groups')
end

% first column of design matrix is a vector of ones
designMatrix = ones(nObs, 1);
%designMatrix = [designMatrix age] %Karin comment
designVars = {'1'};

if opts.groupEffect 
     if groups>1
     for i=1:groups-1
           designMatrix = [designMatrix grouping(:,i)]; % design matrix for model order = 0
           designVars = [designVars; sprintf('grouping_%d', i)];
     end
     else
     warning('only one group, no group effect will nbe taken into account');
     opts.groupEffect=0;
     opts.interEffect=0;
    end
end

if opts.mOrder>0 % if the model order is larger than zero
    for iO = 1:mOrder % add additional colums to the design matrix according to the model order
        designMatrix = [designMatrix age.^iO]; % include age in the model
        if iO==1; designVars = [designVars; 'age']; else; designVars = [designVars; sprintf('age_%d', iO)]; end
        if opts.interEffect % if interaction should be considered  
            designMatrix = [designMatrix (age.^iO).*grouping];
            for i=1:groups-1
               if iO==1; designVars = [designVars; sprintf('age_by_grouping_%d', i)]; else; designVars = [designVars; sprintf('age%d_by_grouping_%d', iO, i);]; end
            end
        end
    end
end

% add covariates if there are any
if ~isempty(cov)
    designMatrix=[designMatrix cov];
    
    for iC=1:size(cov,2)
        designVars=[designVars; ['Covariate ' num2str(iC)]];
    end
end

outModel.input.subjID=subjID;
outModel.input.age=age;
outModel.input.grouping=grouping;
outModel.input.data=data;
outModel.input.cov=cov;

outModel.order=mOrder;
outModel.designMatrix=designMatrix;
outModel.designVars=designVars;

modelFun = @(parameterVector, designMatrix) designMatrix*parameterVector';
initialParameterVectorEstimate = robustfit(designMatrix, data, 'bisquare', 4.685, 'off');

 switch opts.mType
    case 'intercept'
        randParams = 1;
    case 'slope'
        if opts.mOrder>0 % random slope model only if age is included in the model
            if opts.groupEffect
                randParams = [1,groups+1];
            else
                randParams = [1,2];
            end
        else
%           warning('no random slope model for constant model estimation, computing random intercept model');
            randParams = 1;
        end
    case 'glm'
        randParams=[];
  end   
 
[outModel.beta, outModel.randCov, outModel.stats] = nlmefit(designMatrix, data, ...
    subjID, [], modelFun, initialParameterVectorEstimate,'REParamsSelect', randParams,...
    'Options',statset('MaxIter',500),'RefineBeta0',false);
