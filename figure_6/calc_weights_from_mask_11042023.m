function [n_weights_max_abs,n_weights_mean] = calc_weights_from_mask_11042023(Wout,stat)
% 11/08/2022

% Calculate neurons weights from pixels using neuron masks.
% Need to downsample neruon masks

num_neurons = length(stat);
num_x = size(Wout,2);
n_weights_max_abs = zeros(num_neurons,num_x);
n_weights_mean = zeros(num_neurons,num_x);
%% 

for n = 1:num_neurons
    cur_stat = stat{n};
    cur_footprint = zeros(512,512);
    % Downsample the masks
    for j = 1:length(cur_stat.xpix)
        cur_footprint(cur_stat.xpix(j)+1,cur_stat.ypix(j)+1) = 1;
    end
    down_footprint = zeros(128,128);
    for i = 1:128
        for j = 1:128
            if sum(cur_footprint((i-1)*4+1:i*4,(j-1)*4+1:j*4),'all')>0
                down_footprint(i,j) = 1;
            end
        end
    end
    % Caclulate weights in mask
    for i = 1:num_x
        cur_w = reshape(Wout(1:end-1,i),128,128);
        cur_n_w = cur_w(down_footprint==1);
        if ~isempty(cur_n_w)
            n_weights_max_abs(n,i) = max(abs(cur_n_w));
            n_weights_mean(n,i) = mean(cur_n_w);
        end
    end

end