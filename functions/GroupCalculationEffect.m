function results = GroupCalculationEffect(outModelVect)
i=0;
for iM=1:size(outModelVect,1)
if isfield(outModelVect{iM,1}, 'groupEffect')
    if iM==1
    effectSizeGroup=table(iM, outModelVect{1,1}.groupEffect.Chi2, outModelVect{1,1}.groupEffect.dof_diff, outModelVect{1,1}.groupEffect.p);
    effectSizeGroup.Properties.VariableNames={'Number of model','Chi square statistics','Degree of fredoom','p-value'};
    else 
    effectSizeGroupAdd=table(iM, outModelVect{1,1}.groupEffect.Chi2, outModelVect{1,1}.groupEffect.dof_diff, outModelVect{1,1}.groupEffect.p);
    effectSizeGroupAdd.Properties.VariableNames={'Number of model','Chi square statistics','Degree of fredoom','p-value'};
    effectSizeGroup=[effectSizeGroup;effectSizeGroupAdd];   
    end
    results=effectSizeGroup;
    i=i+1;
end
end

if i==0
    results=[];
end


