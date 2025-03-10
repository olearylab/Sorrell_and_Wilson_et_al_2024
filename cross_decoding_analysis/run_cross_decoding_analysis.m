% run cross decoding analysis

% These need to be updated if this script is moved
addpath("../shared_functions")
addpath("../figure_5")

% Set up workspace
virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");
virmen_train_cell = importdata("../saved_files/virmen_train_cell_5.mat");
tbt_train_cell = importdata("../saved_files/tbt_train_cell_5.mat");

nbins = 50;
offsets = [1.4953,1.4961;1.4953,1.4961;1.4953,1.4961;1.4804,1.4844;1.4804,1.4844];

% Calculate decoding accuracy for decoding corrective angular velocities
% using view angle decoder
[all_res_mat] = calc_va_cross_decoding_all(virmen_cell,tbt_cell,virmen_train_cell,tbt_train_cell,nbins,offsets);
