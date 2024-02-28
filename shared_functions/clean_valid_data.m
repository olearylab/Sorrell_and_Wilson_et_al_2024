function [cleaned_valid] = clean_valid_data(ITI)

% function for removing erroneous valid samples from within ITI.
% uses fact that ITI ends when going from +ve to 0, and erroneous samples
% appear when it goes -1 to 0 to 1 (due to binning).

% output is 1s for valid, 0s for invalid

cleaned_valid = ITI == 0;
in_ITI = ~cleaned_valid(1);
for i = 2:length(ITI)
    if ~in_ITI
        cleaned_valid(i) = cleaned_valid(i);
        if cleaned_valid(i) ~= 1
            in_ITI = true;
        end
    else
        if (ITI(i)==0) && ITI(i-1)>0
            in_ITI = false;
            cleaned_valid(i) = 1;
        else
            cleaned_valid(i) = 0;
        end
    end
end