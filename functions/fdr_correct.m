function outModelVect = fdr_correct(outModelVect,alpha)
% fdr_correct - function that takes all p-values in outModelVect, corrects
% them and saves them in the outModelVect

for iM=1:length(outModelVect)
    if isfield(outModelVect{iM},'groupEffect')
        p_gr(iM,1)=outModelVect{iM}.groupEffect.p;
    end
    if isfield(outModelVect{iM},'interEffect')
        p_int(iM,1)=outModelVect{iM}.interEffect.p;
    end
end

if exist('p_gr','var')
    [h_gr_corr, p_gr_crit, ~, p_gr_corr]=fdr_bh(p_gr,alpha);
end
if exist('p_int','var')
    [h_int_corr, p_int_crit, ~, p_int_corr]=fdr_bh(p_int,alpha);
end

for iM=1:length(outModelVect)
    if isfield(outModelVect{iM},'groupEffect')
        outModelVect{iM}.groupEffect.p=p_gr_corr(iM);
        outModelVect{iM}.groupEffect.h=h_gr_corr(iM);
    end
    if isfield(outModelVect{iM},'interEffect')
        outModelVect{iM}.interEffect.p=p_int_corr(iM);
        outModelVect{iM}.interEffect.h=h_int_corr(iM);
    end
end