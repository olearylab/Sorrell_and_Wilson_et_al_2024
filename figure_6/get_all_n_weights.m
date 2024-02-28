function [n_W_cell] = get_all_n_weights(Wout_cell,stat_cell,is_offline,w_ind)
% 27/07/2023

% Get all the weights in neurons, for all sessions

num_mice = size(stat_cell,1);
num_days = size(stat_cell,2);

n_W_cell = cell(num_mice,num_days);

for m = 1:num_mice
    Wout = squeeze(Wout_cell{m}(:,w_ind));
    if is_offline
        if m < 4
            % Transpose needed for offline weights for m1-3
            ww = reshape(Wout(1:end-1),128,128)';
            Wout = ww(:);
            Wout = [Wout;nan];
        end
    end
    for d = 1:num_days
        if ~isempty(stat_cell{m,d})
            % Calculate mean weight within each neuron mask
            [~,n_weights_mean] = calc_weights_from_mask_11042023(Wout,stat_cell{m,d});
            n_W_cell{m,d} = n_weights_mean;
        end
    end
end