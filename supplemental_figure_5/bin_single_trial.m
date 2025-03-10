function [binned_vec] = bin_single_trial(to_bin,binned,nbins)
% 30/09/2024

% function for binning single trial data. 
binned_vec = nan.*ones(nbins,1);
for j = 1:nbins
    if sum(binned==j) == 1
        binned_vec(j) = to_bin(binned==j);
    else
        binned_vec(j) = mean(to_bin(binned==j));
    end
end