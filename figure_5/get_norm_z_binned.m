function [z_binned] = get_norm_z_binned(zdata,virmen_data,tbt_details,model_params,normalise_z,nbins)
% 09/09/2023

% normalise filtered z by ball trials. Then bin and bootstrap

types_vec = [1,4,7,10];
max_types = max(tbt_details(1,:));
types_vec = types_vec(1:max_types);
%% remove invalid data
ITI = virmen_data(8,:);
cleaned_valid = clean_valid_data(ITI);

virmen_data = virmen_data(:,cleaned_valid);
%% Preprocess and bin data

% initialise_params;
if ndims(zdata) == 2
    model_params.spatial = false;
    create_reg = false;
else
    create_reg = true;
end
model_params.reg_images = false;
disp("Pre-Processing")
[zdatafilt, train_mean, Rfixed] = preprocess_images_05_2021(zdata,model_params,[],[],cleaned_valid,create_reg);
disp("Pre-Processing Complete")   
% zdatafilt has invalid data also removed

%% Normalise (z-score) post filtering, before bootstrapping.
% Normalise filtered data with invalid data removed
% normalise using ball trial activity
if normalise_z 
    [zdatafilt] = normalise_z_by_ball(zdatafilt,virmen_data,tbt_details);
end

xnum = 6;
linearise_x = true;
% z_binned is trials x nbins x zdim
[z_binned, centres] = bin_neural_data_general_any(virmen_data,zdatafilt,xnum,linearise_x,nbins);

