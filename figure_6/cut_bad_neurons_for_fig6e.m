% 13/01/2024

% Functioning for removing neurons with invalid data
% 2 neurons were identified, one with many 0s, the other with only integer
% values, neither are trustworthy.

% z_cell should be before cut, to help with indexing the neurons. Should be
% neuron 192 from {3,2} and 200 from {5,3}

% remove from all_full_results_both
all_full_results_cut = all_full_results_both;
temp = all_full_results_both{3,2};
temp_z = z_cell_CNN{3,2};
temp_cut = temp(temp_z(1,:)~=0,:,:);
all_full_results_cut{3,2} = temp_cut;

temp = all_full_results_both{5,3};
temp_z = z_cell_CNN{5,3};
temp_cut = temp(temp_z(1,:)~=162,:,:);
all_full_results_cut{5,3} = temp_cut;

% remove from sig_check_all
sig_results_all_cut = sig_results_all;
temp = sig_results_all{3,2};
temp_z = z_cell_CNN{3,2};
temp_cut = temp(:,temp_z(1,:)~=0);
sig_results_all_cut{3,2} = temp_cut;

temp = sig_results_all{5,3};
temp_z = z_cell_CNN{5,3};
temp_cut = temp(:,temp_z(1,:)~=162);
sig_results_all_cut{5,3} = temp_cut;

% remove from sig_results_all
sig_check_all_cut = sig_check_all;
temp = sig_check_all{3,2};
temp_z = z_cell_CNN{3,2};
temp_cut = temp(temp_z(1,:)~=0);
sig_check_all_cut{3,2} = temp_cut;

temp = sig_check_all{5,3};
temp_z = z_cell_CNN{5,3};
temp_cut = temp(temp_z(1,:)~=162);
sig_check_all_cut{5,3} = temp_cut;

% remove from sig_summary_cell
sig_summary_cell_cut = sig_summary_cell;
temp = sig_summary_cell{3,2};
temp_z = z_cell_CNN{3,2};
temp_cut = temp(:,:,temp_z(1,:)~=0);
sig_summary_cell_cut{3,2} = temp_cut;

temp = sig_summary_cell{5,3};
temp_z = z_cell_CNN{5,3};
temp_cut = temp(:,:,temp_z(1,:)~=162);
sig_summary_cell_cut{5,3} = temp_cut;

% remove from full_all_params_means
full_all_params_means_cut = full_all_params_means;
temp = full_all_params_means{3,2};
temp_z = z_cell_CNN{3,2};
temp_cut = temp(:,temp_z(1,:)~=0,:);
full_all_params_means_cut{3,2} = temp_cut;

temp = full_all_params_means{5,3};
temp_z = z_cell_CNN{5,3};
temp_cut = temp(:,temp_z(1,:)~=162,:);
full_all_params_means_cut{5,3} = temp_cut;
