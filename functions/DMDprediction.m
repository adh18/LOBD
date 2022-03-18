function [pred] = DMDprediction(dmdX, dmdT, cfs)
%DMDPREDICTION Summary of this function goes here
%   Detailed explanation goes here
cfs = cfs(:);
pred = dmdX * (cfs .* dmdT');
end

