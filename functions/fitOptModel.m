function outModelVect = fitOptModel(input, opts)
% fitOptModel - function to fit a random intersect models of order 0 to 3 
% to the input data and return the best model according to the BIC
%
% Syntax:  model = fitOptModel(input, opts)
%
% Inputs:
%    input - structure containing the input information:
%           .subjID - vector of subject IDs (double, #sub x 1)
%           .age - vector of ages (double, #sub x 1)
%           .grouping - vector with grouping information (double/logical)
%           .data - matrix with data to model (double, #sub x #models)
%           [.cov] - matrix with covariates to consider
%    opts - structure with model fitting options
%           [.vertID] - columns of input.data to analyze
%           [.modelNames] - names of models, default = index
%           .orders - vector of model orders to check for
%           [.alpha] - limit for significance check, default = 0.05
%           [.mType] - 'intercept' (default), random intersect model
%                      'slope' for random slope model
%
% Outputs:
%    outModelVect - vector of structs with fitted optimum models (#models x 1)
%           .mName - name of model according to opts
%           .input - estimation inputs
%           .order - order of resulting model
%           .designMatrix - design matrix of estimated model
%           .designVars - names of columns of the design matrix
%           .beta - vector of estimated fixed effects (betas)
%           .randCov - covariance matrix of estimated random effect
%           .stats  - structure with model fitting statistics (see nlmefit
%                   for details)
%           .groupEffect - structure with significance of group difference 
%               (.h and .p) and estimated reduced model (.reducedModel)
%               without grouping variable
%           .interEffect - structure with significance of age*group 
%               interaction effect (.h and .p) and estimated reduced model 
%               (.reducedModel) without interaction
%
% Example: 
%

% Author: Daniela Zoeller
% Medical Image Processing Lab, EPFL/UniGe
% Developmental Imaging and Psychopathology Lab, UniGe
% v1.0 3.8.2016 DZ - initial version based on initial code of Kadir Mutlu

%% set up defaults
if ~isfield(opts,'alpha')
    opts.alpha=0.05;
end

if ~isfield(opts,'vertID')
    opts.vertID=1:size(input.data,2);
end

if ~isfield(opts,'mType')
    opts.mType='intercept';
end

if ~isfield(input,'cov')
    input.cov=[];
else
    input.cov=input.cov-repmat(mean(input.cov),size(input.cov,1),1); % center covariates to avoid confounds
end


%% initialize
nModels=length(opts.vertID);
nOrders=length(opts.orders);

%% estimate models and return best fit according to BIC
for iM = 1:nModels
    fprintf('\nModel %d : ',opts.vertID(iM));
    iV=opts.vertID(iM);
    dataVect=input.data(:,iV);
    
    if any(isnan(dataVect))
        warning('Your data vector is not complete! Running analysis without missing data points');
    end
    dataID=~isnan(dataVect);
    
    if ~any(dataVect)
        if isfield(opts,'modelNames')
            outModelVect{iM,1}.mName = opts.modelNames{iM};
        else
            outModelVect{iM,1}.mName = num2str(iV);
        end
        continue;
    end
    
    for iO = 1:nOrders
        fprintf('%d ',opts.orders(iO));
        estOpts.mOrder=opts.orders(iO);
        estOpts.groupEffect=1;
        estOpts.interEffect=1;
        estOpts.mType=opts.mType;
        % estimate temporary model
        if ~isempty(input.cov)
            tmpModel = estimateModel(input.subjID(dataID,:), input.grouping(dataID,:),...
                input.age(dataID,:), input.cov(dataID,:), dataVect(dataID,:), estOpts);
        else
            tmpModel = estimateModel(input.subjID(dataID,:), input.grouping(dataID,:),...
                input.age(dataID,:), [], dataVect(dataID,:), estOpts);
        end
        if iO==1
            outModelVect{iM,1}=tmpModel; % save constant model
        elseif tmpModel.stats.bic < outModelVect{iM}.stats.bic-2
            outModelVect{iM,1}=tmpModel; % keep higher order model if BIC decreases
        else
            break; % if bic does not decrease sufficiently enough, stop increasing model order (not necessary to compute higher order models)
        end
    end
    
    if isfield(opts,'modelNames')
        outModelVect{iM,1}.mName = opts.modelNames{iM};
    else
        outModelVect{iM,1}.mName = num2str(iV);
    end
    
    
    fprintf('\nFinal model order: %d ',outModelVect{iM}.order);
    
    if length(unique(input.grouping))>1 % if there is more than one group
        %% likelihood ratio test for significant group effect
        estOpts.groupEffect=0;
        estOpts.interEffect=0;
        estOpts.mOrder=outModelVect{iM}.order;
        estOpts.mType=opts.mType;
        if ~isempty(input.cov)
            reducedModel = estimateModel(input.subjID(dataID,:), input.grouping(dataID,:),...
                input.age(dataID,:), input.cov(dataID,:), dataVect(dataID,:), estOpts);
        else
            reducedModel = estimateModel(input.subjID(dataID,:), input.grouping(dataID,:),...
                input.age(dataID,:), [], dataVect(dataID,:), estOpts);
        end

        % degrees of freedom of the chi-squared distribution for the likelihood ratio test is de difference dofs of every model
        dof = length(outModelVect{iM}.beta)-length(reducedModel.beta);
        [outModelVect{iM}.groupEffect.h outModelVect{iM}.groupEffect.p outModelVect{iM}.groupEffect.dof_diff outModelVect{iM}.groupEffect.Chi2] = likelihoodratiotest(outModelVect{iM}.stats.logl, reducedModel.stats.logl, dof, opts.alpha);
        outModelVect{iM}.groupEffect.reducedModel = reducedModel;
    
        %% likelihood ratio test for significant age*group interaction
        if outModelVect{iM}.order
            estOpts.groupEffect=1;
            estOpts.interEffect=0;
            estOpts.mOrder=outModelVect{iM}.order;
            estOpts.mType=opts.mType;
            if ~isempty(input.cov)
                reducedModel = estimateModel(input.subjID(dataID,:), input.grouping(dataID,:),...
                    input.age(dataID,:), input.cov(dataID,:), dataVect(dataID,:), estOpts);
            else
                reducedModel = estimateModel(input.subjID(dataID,:), input.grouping(dataID,:),...
                    input.age(dataID,:), [], dataVect(dataID,:), estOpts);
            end
            dof = length(outModelVect{iM}.beta)-length(reducedModel.beta);
            [outModelVect{iM}.interEffect.h outModelVect{iM}.interEffect.p outModelVect{iM}.interEffect.dof_diff outModelVect{iM}.interEffect.Chi2] = likelihoodratiotest(outModelVect{iM}.stats.logl, reducedModel.stats.logl, dof, opts.alpha);
            outModelVect{iM}.interEffect.reducedModel = reducedModel;
        end
    end
end

end



