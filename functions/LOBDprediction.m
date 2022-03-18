function [pred] = LOBDprediction(lobd, cfs, repeats)
%LOBDPREDICTION Summary of this function goes here
%   Detailed explanation goes here
if nargin == 2
    repeats = 1;
end

cfs = cfs(:);
pred = lobd.factors{1}*diag(cfs)*conj(lobd.factors{2}');

for i = 2:repeats
   tmpcfs = lobd.factors{1}'*pred(:, end);
   %tmpcfs = lobd.factors{1}'*pred(:, end) ./ conj(lobd.factors{2}(1, :)');   % orthogonal projection
   tmppred = lobd.factors{1}*diag(tmpcfs)*conj(lobd.factors{2}');
   pred = cat(2, pred, tmppred(:, 2:end)); 
end

end

