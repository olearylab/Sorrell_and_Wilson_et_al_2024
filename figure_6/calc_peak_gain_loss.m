function [sig_summary_cell] = calc_peak_gain_loss(sig_results_all)
% 24/06/2022

% For each neuron, determine whether it is tuned to both ball and bmi, 
% gained in bmi, lost in bmi. Separately for left and right trials

num_mice = size(sig_results_all,1);
num_days = size(sig_results_all,2);
%%

sig_summary_cell = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(sig_results_all{m,d})
            
            cur_sig = sig_results_all{m,d};
            % [direction, both/gain/loss, neurons]
            sig_summary = zeros(2,3,size(cur_sig,2));
            
            for n = 1:size(cur_sig,2)
                % both left
                if (cur_sig(1,n) == 1) && (cur_sig(3,n) == 1)
                    sig_summary(1,1,n) = 1;
                end
                % both right
                if (cur_sig(2,n) == 1) && (cur_sig(4,n) == 1)
                    sig_summary(2,1,n) = 1;
                end
                % gain left
                if (cur_sig(1,n) == 0) && (cur_sig(3,n) == 1)
                    sig_summary(1,2,n) = 1;
                end
                % lose left
                if (cur_sig(1,n) == 1) && (cur_sig(3,n) == 0)
                    sig_summary(1,3,n) = 1;
                end
                % gain right
                if (cur_sig(2,n) == 0) && (cur_sig(4,n) == 1)
                    sig_summary(2,2,n) = 1;
                end
                % lose right
                if (cur_sig(2,n) == 1) && (cur_sig(4,n) == 0)
                    sig_summary(2,3,n) = 1;
                end
            end
            sig_summary_cell{m,d} = sig_summary;
        end
    end
end