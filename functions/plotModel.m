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


%%
% model=outModelVect{iM};
% plotOpts= plotOpts;
% figurehandle=fig;
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

nCov = plotOpts.nCov;

%% check which model to plot 
% (reduced model if group/interaction effects are not significant)
% define the groups
grouping = model.input.grouping;
if size(grouping,2)==0 
    groups=1;
    data{1}=model.input.data;
    age{1}=model.input.age;
    subjects{1}=model.input.subjID;
    group_logic(:,1)=repmat(1,size(model.designMatrix,1),1);
    else
    groups=size(grouping,2)+1;
end

if groups==2
    data{1}=model.input.data(logical(grouping));
    age{1}=model.input.age(logical(grouping));
    subjects{1}=model.input.subjID(logical(grouping));

    data{2}=model.input.data(logical(~grouping));
    age{2}=model.input.age(logical(~grouping));
    subjects{2}=model.input.subjID(~logical(grouping)); 
    
    group_logic=[logical(grouping), logical(~grouping)];
end

if groups==3
    for iG = 1:groups-1
    data{iG}=model.input.data(logical(grouping(:,iG)));
    age{iG}=model.input.age(logical(grouping(:,iG)));
    subjects{iG}=model.input.subjID(logical(grouping(:,iG)));
    end
    for i=1:size(grouping,1)
       if grouping(i,1)==0 & grouping(i,2)==0
       idx_3group(i)=1;
       else
       idx_3group(i)=0;
       end
    end
    for iG=groups
    data{iG}=model.input.data(logical(idx_3group));
    age{iG}=model.input.age(logical(idx_3group));
    subjects{iG}=model.input.subjID(logical(idx_3group)); 
    end

    group_logic=[logical(grouping(:,1)), logical(grouping(:,2)), logical(idx_3group')];
end

figure(fh);
if strcmp(plotOpts.plotType,'redGrp') && (~isfield(model,'groupEffect') || ~model.groupEffect.h) % no significant group effect
    plModel = model.groupEffect.reducedModel;
    
    % get model parameters
    params = plModel.beta(1:end-nCov)';
    
    % plot data
    for iG = 1:groups
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
    for iG=1:groups
        pl(iG) = plotDataAndFit(age{iG},data{iG},subjects{iG},[],plotCol{iG});
    end
    
    plModel = model.interEffect.reducedModel;
    
    for iG = 1:groups
        R{iG} = plModel.stats.iwres(group_logic(:,iG));
    end
 
    if groups==2
        if model.order==1
             params{1}(1) = plModel.beta(1)+plModel.beta(2);
             params{1}(2) = plModel.beta(3);
             params{2}(1) = plModel.beta(1);
             params{2}(2) = plModel.beta(3);
        end
        if  model.order==2
             params{1}(1) = plModel.beta(1)+plModel.beta(2);
             params{1}(2) = plModel.beta(3);
             params{1}(3) = plModel.beta(4);
             params{2}(1) = plModel.beta(1);
             params{2}(2) = plModel.beta(3);
             params{2}(3) = plModel.beta(4);   
        end
        if model.order==3
             params{1}(1) = plModel.beta(1)+plModel.beta(2);
             params{1}(2) = plModel.beta(3);
             params{1}(3) = plModel.beta(4);
             params{1}(4) = plModel.beta(5);
             params{2}(1) = plModel.beta(1);
             params{2}(2) = plModel.beta(3);
             params{2}(3) = plModel.beta(4);   
             params{2}(4) = plModel.beta(5);        
        end
    end
    if groups==3
        if model.order==1
             params{1}(1) = plModel.beta(1)+plModel.beta(2);
             params{2}(1) = plModel.beta(1)+plModel.beta(3);
             params{3}(1) = plModel.beta(1);
             params{1}(2) = plModel.beta(4);
             params{2}(2) = plModel.beta(4);
             params{3}(2) = plModel.beta(4);
        end
        if model.order==2
             params{1}(1) = plModel.beta(1)+plModel.beta(2);
             params{2}(1) = plModel.beta(1)+plModel.beta(3);
             params{3}(1) = plModel.beta(1);
             params{1}(2) = plModel.beta(4);
             params{2}(2) = plModel.beta(4);
             params{3}(2) = plModel.beta(4);
             params{1}(3) = plModel.beta(5);
             params{2}(3) = plModel.beta(5);
             params{3}(3) = plModel.beta(5);
             
        end
        if model.order==3
             params{1}(1) = plModel.beta(1)+plModel.beta(2);
             params{2}(1) = plModel.beta(1)+plModel.beta(3);
             params{3}(1) = plModel.beta(1);
             params{1}(2) = plModel.beta(4);
             params{2}(2) = plModel.beta(4);
             params{3}(2) = plModel.beta(4);
             params{1}(3) = plModel.beta(5);
             params{2}(3) = plModel.beta(5);
             params{3}(3) = plModel.beta(5);
             params{1}(4) = plModel.beta(6);
             params{2}(4) = plModel.beta(6);
             params{3}(4) = plModel.beta(6);
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
        for iG=1:groups
            var_id=group_logic(:,iG);
            y{iG}=ypred(var_id(sortID));
            delt{iG}=delta(var_id(sortID));
            ageSort{iG}=ageS(var_id(sortID));
            
            plotDataAndFit(age{iG},[],[],params{iG},plotCol{iG},ageSort{iG},y{iG},delt{iG});
        end
    else
        % plot fitted models on top
        for iG=1:groups         
            plotDataAndFit(age{iG},[],[],params{iG},plotCol{iG});
        end
    end
    
    plottedModel = 'noInteraction';
else
    % plot data
    for iG=1:groups
        pl(iG) = plotDataAndFit(age{iG},data{iG},subjects{iG},[],plotCol{iG});
    end
    
    plModel = model;  
    if groups>1
      for iG = 1:groups
            R{iG} = plModel.stats.iwres(group_logic(:,iG));
      end
      
      if plModel.order==0 
          if groups==2 
          params{1}(1)=plModel.beta(1)+plModel.beta(2);
          params{2}(1)=plModel.beta(1);
          end
          if groups==3
          params{1}(1)=plModel.beta(1)+plModel.beta(2);
          params{2}(1)=plModel.beta(1)+plModel.beta(3);
          params{3}(1)=plModel.beta(1);
          end
      end
      
      if plModel.order==1
          if groups==2 
          params{1}(1)=plModel.beta(1)+plModel.beta(2);
          params{1}(2)=plModel.beta(3)+plModel.beta(4);
          params{2}(1)=plModel.beta(1);
          params{2}(2)=plModel.beta(3);
          end
          if groups==3
          params{1}(1)=plModel.beta(1)+plModel.beta(2);
          params{1}(2)=plModel.beta(4)+plModel.beta(5);
          params{2}(1)=plModel.beta(1)+plModel.beta(3);
          params{2}(2)=plModel.beta(4)+plModel.beta(6);
          params{3}(1)=plModel.beta(1);
          params{3}(2)=plModel.beta(4);
          end
      end
      
      if plModel.order==2 
          if groups==2
          params{1}(1)=plModel.beta(1)+plModel.beta(2);
          params{1}(2)=plModel.beta(3)+plModel.beta(4);
          params{1}(3)=plModel.beta(5)+plModel.beta(6);
          params{2}(1)=plModel.beta(1);
          params{2}(2)=plModel.beta(3);
          params{2}(3)=plModel.beta(5);
          end
          if groups==3
          params{1}(1)=plModel.beta(1)+plModel.beta(2);
          params{1}(2)=plModel.beta(4)+plModel.beta(5);
          params{1}(3)=plModel.beta(7)+plModel.beta(8);
          params{2}(1)=plModel.beta(1)+plModel.beta(3);
          params{2}(2)=plModel.beta(4)+plModel.beta(6);
          params{2}(3)=plModel.beta(7)+plModel.beta(9);
          params{3}(1)=plModel.beta(1);
          params{3}(2)=plModel.beta(4);
          params{3}(3)=plModel.beta(7);
          end
      end
      
      if plModel.order==3
         if groups==2
         params{1}(1)=plModel.beta(1)+plModel.beta(2);
         params{1}(2)=plModel.beta(3)+plModel.beta(4);
         params{1}(3)=plModel.beta(5)+plModel.beta(6);
         params{1}(4)=plModel.beta(7)+plModel.beta(8);
         params{2}(1)=plModel.beta(1);
         params{2}(2)=plModel.beta(3);
         params{2}(3)=plModel.beta(5);
         params{2}(4)=plModel.beta(7);  
         end
         if groups==3
         params{1}(1)=plModel.beta(1)+plModel.beta(2);
         params{1}(2)=plModel.beta(4)+plModel.beta(5);
         params{1}(3)=plModel.beta(7)+plModel.beta(8);
         params{1}(4)=plModel.beta(10)+plModel.beta(11);
         params{2}(1)=plModel.beta(1)+plModel.beta(3);
         params{2}(2)=plModel.beta(4)+plModel.beta(6);
         params{2}(3)=plModel.beta(7)+plModel.beta(9);
         params{2}(4)=plModel.beta(10)+plModel.beta(12);
         params{3}(1)=plModel.beta(1);
         params{3}(2)=plModel.beta(4);
         params{3}(3)=plModel.beta(7);
         params{3}(4)=plModel.beta(10);
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
        for iG=1:groups
            var_id=group_logic(:,iG);
            y{iG}=ypred(var_id(sortID));
            delt{iG}=delta(var_id(sortID));
            ageSort{iG}=ageS(var_id(sortID));
            
            plotDataAndFit(age{iG},[],[],params{iG},plotCol{iG},ageSort{iG},y{iG},delt{iG});
        end
    else
        % plot fitted models on top
        for iG=1:groups         
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
    for iG=1:groups
        legTxt{iG}=['Group ' num2str(iG)];
    end
    legend(pl,legTxt);
end
