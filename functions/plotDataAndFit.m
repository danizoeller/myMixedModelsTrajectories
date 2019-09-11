% function pl = plotDataAndFit(age,data,subjects,params,plotColor,R,CovB,desMat)
function pl = plotDataAndFit(age,data,subjects,params,plotColor,ageS,ypred,delta)
% plotDataAndFit - function plot longitudinal data points and a fitted
% model curve
%
% Syntax:  plotDataAndFit(age,data,subjects,params,plotSpec)
%
% Inputs:
%    age - observation time points (to plot on x axis)
%    data - model input data (to plot on y axis), if empty no data is plotted
%    subjects - vector of subject IDs (double, #sub x 1)
%    params - model parameters, if empty no curve is plotted
%    plotColor - plotting color
%
% Outputs:
%    graph plotted on current figure
%
% Example: 
%

% Author: Daniela Zoeller
% Medical Image Processing Lab, EPFL/UniGe
% Developmental Imaging and Psychopathology Lab, UniGe
% v1.0 3.8.2016 DZ - initial version based on initial code of Kadir Mutlu

hold on;

% if size(params,1)>size(params,2)
%     warning('parameters have to be given in a row vector, flipping dimensions');
%     params=params';
% end

%% plot data
if ~isempty(data)
    pl = plot(age,data,'.','markersize',21,'Color',plotColor);pl.Color(4)=0.5; 
    
    % plot lines per subject
    uniqueSub=unique(subjects);
    nSub=length(uniqueSub);
    for iS = 1:nSub
        sID=uniqueSub(iS);
        pt=plot(age(subjects==sID),data(subjects==sID),'-.','Color',plotColor);pt.Color(4)=0.5;
    end
end

%% plot fitted model
if ~isempty(params)
    ageVec=min(age):0.1:max(age);
    plot(ageVec,polyval(flip(params),ageVec),'LineWidth',3,'Color',plotColor*0.8);
    
%     if nargin == 8
%         % confidecte intervals
%         modelfun = @(parameterVector, designMatrix) designMatrix*parameterVector';
%         [ageS,sortID]=sort(age);
%         beta=params;
%         desMat=desMat(sortID,:);
%         R=R(sortID);
%         [ypred,delta] = nlpredci(modelfun,desMat,beta,R,'Covar',CovB);
%         % confplot(ageS,ypred,delta);
%         shadedErrorBar(ageS,ypred,delta,plotColor,0.5)
%     end

    if nargin == 8 && exist('ypred','var') && ~isempty(ypred)
        shadedErrorBar(ageS,ypred,delta,{'Color',plotColor*0.8},0.5)
    end
end
