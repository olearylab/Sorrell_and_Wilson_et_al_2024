function [final_footprints] = downsample_neuron_masks(stat)
% 06/10/2023

% Downsample neuron masks from 512 x 512 to 128 x 128

% If any pixel in each 4x4 block is in full neuron mask, then assign
% donwsampled pixel to the mask

num_neurons = length(stat);
final_footprints = zeros(128,128);
for n = 1:num_neurons
    cur_stat = stat{n};
    cur_footprint = zeros(512,512);
    for j = 1:length(cur_stat.xpix)
        cur_footprint(cur_stat.xpix(j)+1,cur_stat.ypix(j)+1) = 1;
    end
    down_footprint = zeros(128,128);
    for i = 1:128
        for j = 1:128
            if sum(cur_footprint((i-1)*4+1:i*4,(j-1)*4+1:j*4),'all')>0
                down_footprint(i,j) = 1;
                final_footprints(i,j) = 1;
            end
        end
    end
     
end