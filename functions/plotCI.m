function plotCI(age,params,R,CovB,desMat)

hold on

%% confidence intervals
% confidecte intervals
modelfun = @(parameterVector, designMatrix) designMatrix*parameterVector';
[ageS,sortID]=sort(age);
beta=params;
desMat=desMat(sortID,:);
R=R(sortID);
[ypred,delta] = nlpredci(modelfun,desMat,beta,R,'Covar',CovB);
% confplot(ageS,ypred,delta);
shadedErrorBar(ageS,ypred,delta,[],0.5)
