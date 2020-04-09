function results = InterCalculationEffect(outModelVect)
i=0;
for iM=1:size(outModelVect,1)
if isfield(outModelVect{iM,1}, 'interEffect')
    if iM==1
    effectSizeInter=table(iM, outModelVect{1,1}.interEffect.Chi2, outModelVect{1,1}.interEffect.dof_diff, outModelVect{1,1}.interEffect.p);
    effectSizeInter.Properties.VariableNames={'Number of model','Chi square statistics','Degree of fredoom','p-value'};
    else 
    effectSizeInterAdd=table(iM, outModelVect{1,1}.interEffect.Chi2, outModelVect{1,1}.interEffect.dof_diff, outModelVect{1,1}.interEffect.p);
    effectSizeInterAdd.Properties.VariableNames={'Number of model','Chi square statistics','Degree of fredoom','p-value'};
    effectSizeInter=[effectSizeInter;effectSizeInterAdd];
    end
    results=effectSizeInter;
    i=i+1;
end
end

if i==0
    results=[];
end
