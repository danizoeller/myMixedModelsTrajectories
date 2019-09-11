function sig = com_sig(p_val)
% PURPOSE: To compute the significance
%---------------------------------------------------
% USAGE: sig = com_sig(p_val)
% where:
% p_val = P-values
%---------------------------------------------------
% RETURNS: a vector of output arguments composed of:
% sig = Significances
% --------------------------------------------------
% SEE ALSO: f(results)
%---------------------------------------------------
% REFERENCES: 
%---------------------------------------------------
% REMARKS: 
%---------------------------------------------------

% Written by: Kadir Mutlu

zeroLogicals = p_val == 0;
p_val(zeroLogicals) = realmin;
sig = -log10(p_val);