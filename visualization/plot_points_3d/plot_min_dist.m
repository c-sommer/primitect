function plot_min_dist(in_folder_pc, in_folder_res, fileprefix)
% PLOT_MIN_DIST plots the distance from points in a point cloud to a given set of geometric primitives
%
% PLOT_MIN_DIST(in_folder_pc, in_folder_red, fileprefix)
%
% Input:
% - in_folder_pc: string with folder containing the input ply file
% - in_folder_res: string with folder containing the result files
% - fileprefix: filename of point cloud without '.ply'
%
% No Output.
%
% See also: MIN_DIST.

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer, Yumin Sun.

plyfile = [in_folder_pc fileprefix '.ply'];
pla_file = [in_folder_res fileprefix '_planes.txt'];
sph_file = [in_folder_res fileprefix '_spheres.txt'];
cyl_file = [in_folder_res fileprefix '_cylinders.txt'];
con_file = [in_folder_res fileprefix '_cones.txt'];

pc = pcread(plyfile);

plas = load(pla_file); % p x 5
sphs = load(sph_file); % s x 5
cyls = load(cyl_file); % cy x 8
cons = load(con_file); % co x 8

planes = plas(:,2:end);
spheres = sphs(:,2:end);
cylinders = cyls(:,2:end);
cones = cons(:,2:end);

eval = min_dist(pc,planes,spheres,cylinders,cones);

pcshow(pc.Location, eval)

gray = .75 * [1 1 1];
set(gcf,'color', gray)
set(gca,'color', gray)
axis off