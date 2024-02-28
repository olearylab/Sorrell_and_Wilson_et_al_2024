function [error_mat,x_binned] = BMI_correct_errors_norm(virmen_data,tbt_details,mean_binned,std_binned,nbins,offsets)
% 09/09/2023

% Function for assessing how often the BMI output is error correcting. Just
% like checking if behaviour is error correcting.

%%
% x_binned is trials x bins x xvec [6,7,15,17];
% Calculate heading deviations
[error_mat,x_binned] = check_error_correction_norm(virmen_data,tbt_details,mean_binned,std_binned,nbins,offsets,false);

% Convert ball voltages into velocities
circum = 64;
V = 0.32;
alpha = -50/75*circum/V;
beta = 0.05*circum/V/2.5;
virmen_data(13,:) = alpha*(virmen_data(13,:)-offsets(1));
virmen_data(15,:) = -beta*(virmen_data(15,:)-offsets(2));

dt = 1/30;
% Differentiate decoded view angle
vav = [virmen_data(17,2:end) - virmen_data(17,1:end-1),0]./dt;
virmen_data = [virmen_data;vav];
vav_ind = size(virmen_data,1);

% Bin behaviour including angular velocity from decoded view angle
linearise_x = true;
x_vec = [6,7,15,17,vav_ind];
virmen_data(7,:) = wrapToPi(virmen_data(7,:));
[x_binned,centres] = bin_kin_data(virmen_data,x_vec,linearise_x,nbins);
