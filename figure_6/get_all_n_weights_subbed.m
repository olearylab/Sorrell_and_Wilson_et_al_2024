function [n_W_cell] = get_all_n_weights_subbed(Wout_cell,stat_cell,w_ind)
% 13/09/2023

% Get all the weights in neurons, 
% For when there are many subsamples and cross validations
% This is very slow to run

num_mice = size(stat_cell,1);
num_days = size(stat_cell,2);

n_W_cell = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(stat_cell{m,d})
            disp("running mouse " + m + " day " + d)
            Wout_all = squeeze(Wout_cell{m,d});
            cur_w = nan(size(Wout_all,1),size(Wout_all,2),length(stat_cell{m,d}));
            
            for s = 1:size(cur_w,1)
                for c = 1:size(cur_w,2)
                    Wout = squeeze(Wout_all(s,c,:,w_ind));
                    % Need transpose for first 3
                    if m < 4
                        ww = reshape(Wout(1:end-1),128,128)';
                        Wout = ww(:);
                        Wout = [Wout;nan];
                    end
                    % Calculate mean weights in neuron masks
                    [~,n_weights_mean] = calc_weights_from_mask_11042023(Wout,stat_cell{m,d});
                    cur_w(s,c,:) = n_weights_mean;
                end
            end
            n_W_cell{m,d} = cur_w;
        end
    end
end
save n_W_cell_corrective_subs.mat n_W_cell -v7.3