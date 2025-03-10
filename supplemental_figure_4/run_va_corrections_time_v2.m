%% Script for plotting view angle corrections - supplement figure 4

% These need to be updated if this script is moved
addpath("../shared_functions")

% Set up workspace
virmen_cell = importdata("../saved_files/virmen_cell_5.mat");
tbt_cell = importdata("../saved_files/tbt_cell_5.mat");

% Plotting over time (rather than binned by position)
% Mouse 1
thresh = 0.25;
ex_md = [1,2];
t_ind = 1;
plot_va_corrections_time_v2(virmen_cell, tbt_cell, ex_md, thresh,t_ind);

% Mouse 2
thresh = 0.25;
ex_md = [2,2];
t_ind = 1;
plot_va_corrections_time_v2(virmen_cell, tbt_cell, ex_md, thresh,t_ind);

% Mouse 3
thresh = 0.25;
ex_md = [3,1];
t_ind = 2;
plot_va_corrections_time_v2(virmen_cell, tbt_cell, ex_md, thresh,t_ind);

% Mouse 4
thresh = 0.25;
ex_md = [4,1];
t_ind = 4;
plot_va_corrections_time_v2(virmen_cell, tbt_cell, ex_md, thresh,t_ind);

% Mouse 5
thresh = 0.25;
ex_md = [5,1];
t_ind = 1;
plot_va_corrections_time_v2(virmen_cell, tbt_cell, ex_md, thresh,t_ind);

% If plotting binned use this instead for each mouse and t_nums as above
% linearise_x = true;
% nbins = 50;
% plot_va_corrections(virmen_cell, tbt_cell, ex_md, linearise_x, nbins, t_nums)