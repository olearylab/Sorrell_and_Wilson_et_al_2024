function [sorted_means,I] = sort_means_04_2022(means,normalise,II)
% If providing an order, sort according to this, otherwise sort according
% to the provided means
if isempty(II)
    [maxval,ind] = max(means);
    [b,I] = sort(ind);
else
    I = II;
end
sorted_means = means(:,I);
% normalise between 0 and 1
if normalise
    for i = 1:size(sorted_means,2)
        sorted_means(:,i) = (sorted_means(:,i)-min(sorted_means(:,i)))/(max(sorted_means(:,i))-min(sorted_means(:,i)));
    end   
end
sorted_means = sorted_means';