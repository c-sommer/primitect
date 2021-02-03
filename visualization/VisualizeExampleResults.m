%% Visualize example results for PrimiTect data

in_folder_pc = '../data/primitect/';
in_folder_res = '../ex_results/primitect/';

prefix = '717416';

% display distance to closest geometric primitive instance
figure
plot_min_dist(in_folder_pc, in_folder_res, prefix)

% display labels of found geometric primitives
figure
plot_labels(in_folder_pc, in_folder_res, prefix)

%% Visualize example results for Redwood data

in_folder_pc = '../data/redwood_lod/cylinders/cut2m/';
in_folder_im = '../data/redwood_lod/cylinders/';
in_folder_res = '../ex_results/redwood_lod/cylinders/';

prefix = '01220';

figure
object_overlay(in_folder_pc, in_folder_res, in_folder_im, prefix, 'cylinder')