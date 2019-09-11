function plotModelsAndSaveResults(outModelVect,plotOpts,saveResults,outDir)
% plotModel - function to plot the fitted mixed effect model
%
% Syntax:  plotModelsAndSaveResults(outModelVect,figLabels,[saveResults],[outDir])
%
% Inputs:
%    outModelVect - cell vector of estimated mixed models
%    plotOpts - struct with plotting color and the axis label and legend 
%       information, if empty, the variable names will be used as labels
%       [.plotType] = 'full' - always plot full model, even if interation or 
%                       intercept are not significantly different
%                'redInter' - plot reduced model without interaction if 
%                       interaction is not significant
%                'redGrp' - plot reduced model without group effect if 
%                       group effect is not significant
%       [.plotting] - 0 for no plot saving (default = 1)
%       [.plotCI] - 0 or 1 (default) indicating if confidence bands should be plotted
%       [.plotCol] - cellstring with plotting colors
%    saveResults - 0: no results will be saved
%                  1: results will be saved in outDir
%                  2: results and figures will be saved in outDir
%    outDir - directory for results saving
%
% Outputs:
%    figure plotting and saving of the results if asked for
%
% Example: 
%

% Author: Daniela Zoeller
% Medical Image Processing Lab, EPFL/UniGe
% Developmental Imaging and Psychopathology Lab, UniGe
% v1.0 3.8.2016 DZ - initial version based on initial code of Kadir Mutlu


if nargin < 3
    saveResults=0;
    outDir=[];
end

if nargin < 4 && saveResults
    outDir = pwd;
    warning('no output directory specified, saving to current directory');
end

if isempty(plotOpts) || ~isfield(plotOpts,'plotType')
    plotOpts.plotType='redGrp';
end

if ~isfield(plotOpts,'plotting')
    plotOpts.plotting=1;
end

if ~isfield(plotOpts,'plotCI')
    plotOpts.plotCI=1;
end

if ~isfield(plotOpts,'figPosition')
    plotOpts.figPosition=[440   488   525   310];
end

%% save results (mat file and table)
nMod=length(outModelVect);
if saveResults
    mkdir(outDir); 
    
    for iM=1:nMod
        % define rows to save in the table
        tab.modelName{iM,1} = outModelVect{iM}.mName;
        
        if isfield(outModelVect{iM},'order')
            tab.modelOrder(iM,1) = outModelVect{iM}.order;
            
            if isfield(outModelVect{iM},'groupEffect')
                tab.pValGroup(iM,1)=outModelVect{iM}.groupEffect.p;
            end
            
            if outModelVect{iM}.order & isfield(outModelVect{iM},'interEffect')
                tab.pValInteraction(iM,1)=outModelVect{iM}.interEffect.p;
            end
            
%             if ~outModelVect{iM}.groupEffect.h % no significant group effect
%                 outModel = outModelVect{iM}.groupEffect.reducedModel;
%             elseif outModelVect{iM}.order && ~outModelVect{iM}.interEffect.h % no significant age by group effect
%                 outModel = outModelVect{iM}.interEffect.reducedModel;
%             else
%                 outModel = outModelVect{iM};
%             end
% 
%             for iB = 1:length(outModel.designVars)
%                 fieldname=genvarname(['beta_' outModel.designVars{iB}]);
%                 tab.(fieldname)(iM,1) = outModel.beta(iB);
%             end
            
            for iB = 1:length(outModelVect{iM}.designVars)
                fieldname=genvarname(['full_beta_' outModelVect{iM}.designVars{iB}]);
                tab.(fieldname)(iM,1) = outModelVect{iM}.beta(iB);
            end
            
            if isfield(outModelVect{iM},'groupEffect')
                for iB = 1:length(outModelVect{iM}.groupEffect.reducedModel.designVars)
                    fieldname=genvarname(['noGroup_beta_' outModelVect{iM}.groupEffect.reducedModel.designVars{iB}]);
                    tab.(fieldname)(iM,1) = outModelVect{iM}.groupEffect.reducedModel.beta(iB);
                end
            end
            
            if outModelVect{iM}.order & isfield(outModelVect{iM},'interEffect')
                for iB = 1:length(outModelVect{iM}.interEffect.reducedModel.designVars)
                    fieldname=genvarname(['noInteraction_beta_' outModelVect{iM}.interEffect.reducedModel.designVars{iB}]);
                    tab.(fieldname)(iM,1) = outModelVect{iM}.interEffect.reducedModel.beta(iB);
                end
            end
        end
    end
    
    % fill fields with zeros
    fieldList = fieldnames(tab);
    for iF = 1:length(fieldList)
        if ~iscell(tab.(fieldList{iF}))
            tab.(fieldList{iF,1})(end+1:nMod,1)=0;
        end
    end
    
    % create result table
    tab = struct2table(tab);

    % output file
    if length(outModelVect)>1
        resFile=fullfile(outDir,['resultTable' outModelVect{1}.mName 'to' outModelVect{end}.mName]);
    else
        resFile=fullfile(outDir,['resultTable' outModelVect{1}.mName]);
    end
    writetable(tab,[resFile '.txt'],'Delimiter','\t');
    
    save([resFile '.mat'],'outModelVect');
    

end


%% loop over all models to plot and save
if plotOpts.plotting
    for iM=1:length(outModelVect)
        if isfield(outModelVect{iM},'order') % if there is a model       
            % plotting model
            fig=figure('Position',plotOpts.figPosition);
            [fig,plottedModel]=plotModel(outModelVect{iM},plotOpts,fig);
            
            [fig2,res]=plotResiduals(outModelVect{iM});
            
            % save figure
            if saveResults == 2
                print(fig,fullfile(outDir,[outModelVect{iM}.mName '_' plottedModel]),'-depsc2');
%                 saveas(fig,fullfile(outDir,[outModelVect{iM}.mName '_' plottedModel]));
                
                print(fig2,fullfile(outDir,[outModelVect{iM}.mName '_' plottedModel '_ResNormplot']),'-depsc2');
%                 saveas(fig2,fullfile(outDir,[outModelVect{iM}.mName '_' plottedModel '_ResNormplot']));
            end
        end
    end
end




