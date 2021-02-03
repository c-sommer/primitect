function plot_labels(in_folder_pc, in_folder_res, fileprefix)
% PLOT_LABELS computes and plots the labels of a point cloud consisteing of different geometric primitives
%
% PLOT_LABELS(in_folder_pc, in_folder_red, fileprefix)
%
% Input:
% - in_folder_pc: string with folder containing the input ply file
% - in_folder_res: string with folder containing the result files
% - fileprefix: filename of point cloud without '.ply'
%
% No Output.

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

labels = generate_labels(in_folder_pc, in_folder_res, fileprefix);

pc = pcread([in_folder_pc, fileprefix, '.ply']);
pcshow_labels(pc, labels)