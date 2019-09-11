function [fh,plottedModel] = plotModel(model,plotOpts,figurehandle)
% plotModel - function to plot the fitted mixed effect model
%
% Syntax:  fh = plotModel(model,[figurehandle])
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
%       [.plotCI] - 0 or 1 (default) indicating if confidence bands should be plotted
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
if nargin < 3 || isempty(figurehandle)
    fh=figure;
else
    fh = figurehandle;
end

if nargin < 2
    plotOpts = struct;
end

if ~isfield(plotOpts,'plotCol')
    plotCol = {[0 0 1],[1 0 0],[0 0 0]};
else
    plotCol = plotOpts.plotCol;
end

if ~isfield(plotOpts,'plotType')
    plotOpts.plotType = 'redGrp';
end

nCov = plotOpts.nCov;%length(model.beta)-2*(model.order+1);

%% check which model to plot 
% (reduced model if group/interaction effects are not significant)
grouping = model.input.grouping;
groups=unique(grouping);

for iG = 1:length(groups)
    data{iG}=model.input.data(grouping==groups(iG));
    age{iG}=model.input.age(grouping==groups(iG));
    subjects{iG}=model.input.subjID(grouping==groups(iG));
end


figure(fh);
if strcmp(plotOpts.plotType,'redGrp') && (~isfield(model,'groupEffect') || ~model.groupEffect.h) % no significant group effect
    plModel = model.groupEffect.reducedModel;
    
    % get model parameters
    params = plModel.beta(1:end-nCov)';
    
    % plot data
    for iG = 1:length(groups)
        pl(iG) = plotDataAndFit(age{iG},data{iG},subjects{iG},[],plotCol{iG});
    end
    
    % plotting fitted model
    if plotOpts.plotCI
        % confidecte intervals preparation
        modelfun = @(parameterVector, designMatrix) designMatrix*parameterVector';
        [ageS,sortID]=sort(model.input.age);
        beta=plModel.beta(1:end-nCov)';
        desMat=plModel.designMatrix(sortID,1:end-nCov);
        R=plModel.stats.iwres(sortID);
        CovB=plModel.stats.covb(1:end-nCov,1:end-nCov);
        [ypred,delta] = nlpredci(modelfun,desMat,beta,R,'Covar',CovB);
        
        plotDataAndFit(model.input.age,[],[],params,plotCol{end},ageS,ypred,delta);
    else
        plotDataAndFit(model.input.age,[],[],params,plotCol{end});
    end

    plottedModel = 'noGroup';
    
elseif model.order && ~strcmp(plotOpts.plotType,'full') && isfield(model,'interEffect') && ~model.interEffect.h % no significant age by group effect
    % plot data
    for iG=1:length(groups)
        pl(iG) = plotDataAndFit(age{iG},data{iG},subjects{iG},[],plotCol{iG});
    end
    
    
    plModel = model.interEffect.reducedModel;
    for iG = 1:length(groups)
        R{iG} = plModel.stats.iwres(grouping==groups(iG));
        % get model parameters
        params{iG}(1)= plModel.beta(1) + plModel.beta(2)*groups(iG);
        for iP = 2:length(plModel.beta)-nCov-1
            params{iG}(iP) = plModel.beta(iP+1);
        end
    end
    
    % plotting fitted model per group
    if plotOpts.plotCI
        % confidecte intervals preparation
        modelfun = @(parameterVector, designMatrix) designMatrix*parameterVector';
        [ageS,sortID]=sort(model.input.age);
        beta=plModel.beta(1:end-nCov)';
        desMat=plModel.designMatrix(sortID,1:end-nCov);
        R=plModel.stats.iwres(sortID);
        CovB=plModel.stats.covb(1:end-nCov,1:end-nCov);
        [ypred,delta] = nlpredci(modelfun,desMat,beta,R,'Covar',CovB);
        
        % plot fitted models with confidence intervals
        for iG=1:length(groups)
            y{iG}=ypred(grouping(sortID)==groups(iG));
            delt{iG}=delta(grouping(sortID)==groups(iG));
            ageSort{iG}=ageS(grouping(sortID)==groups(iG));
            
            plotDataAndFit(age{iG},[],[],params{iG},plotCol{iG},ageSort{iG},y{iG},delt{iG});
        end
    else
        % plot fitted models on top
        for iG=1:length(groups)         
            plotDataAndFit(age{iG},[],[],params{iG},plotCol{iG});
        end
    end
    
    plottedModel = 'noInteraction';
else
    % plot data
    for iG=1:length(groups)
        pl(iG) = plotDataAndFit(age{iG},data{iG},subjects{iG},[],plotCol{iG});
    end
    
    
    
    plModel = model;
    if length(groups)>1
        for iG = 1:length(groups)
            R{iG} = plModel.stats.iwres(grouping==groups(iG));
            % get model parameters
            for iP = 1:plModel.order+1
                params{iG}(iP) = plModel.beta(2*iP-1)+plModel.beta(2*iP)*groups(iG);
            end
        end
    else
        R{1} = plModel.stats.iwres;
        % get model parameters
        params{1} = plModel.beta(1:end-nCov);
    end
    
    % plotting fitted model
    if plotOpts.plotCI
        % confidecte intervals preparation
        modelfun = @(parameterVector, designMatrix) designMatrix*parameterVector';
        [ageS,sortID]=sort(model.input.age);
        beta=plModel.beta(1:end-nCov)';
        desMat=plModel.designMatrix(sortID,1:end-nCov);
        R=plModel.stats.iwres(sortID);
        CovB=plModel.stats.covb(1:end-nCov,1:end-nCov);
        [ypred,delta] = nlpredci(modelfun,desMat,beta,R,'Covar',CovB);
        
        % plot fitted models with confidence intervals
        for iG=1:length(groups)
            y{iG}=ypred(grouping(sortID)==groups(iG));
            delt{iG}=delta(grouping(sortID)==groups(iG));
            ageSort{iG}=ageS(grouping(sortID)==groups(iG));
            
            plotDataAndFit(age{iG},[],[],params{iG},plotCol{iG},ageSort{iG},y{iG},delt{iG});
        end
    else
        % plot fitted models on top
        for iG=1:length(groups)         
            plotDataAndFit(age{iG},[],[],params{iG},plotCol{iG});
        end
    end
    
    plottedModel = 'full';
end

%% plot labels
% title
if isfield(plotOpts,'titletext')
    title(plotOpts.titletext)
else
    if ~isfield(model,'groupEffect')
        title(['model order: ' num2str(model.order)]);
    elseif model.order
        title(['model order: ' num2str(model.order) ...
            ',   p-val group effect: ' num2str(model.groupEffect.p,'%.4f') ...
            ',   p-val interaction: ' num2str(model.interEffect.p,'%.4f')]);
    else
        title(['model order: ' num2str(model.order) ...
            ',   p-val group effect: ' num2str(model.groupEffect.p,'%.4f')]);
    end
end

% xlabel
if isfield(plotOpts,'xLabel')
    xlabel(plotOpts.xLabel);
else
    xlabel('age')
end

% ylabel
if isfield(plotOpts,'yLabel')
    ylabel(plotOpts.yLabel);
else
    ylabel('data');
end

% legend
if isfield(plotOpts,'legTxt')
    legend(pl,plotOpts.legTxt)
else
    for iG=1:length(groups)
        legTxt{iG}=['Group ' num2str(iG)];
    end
    legend(pl,legTxt);
end





