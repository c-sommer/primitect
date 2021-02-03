function [labels, types] = generate_labels(input_folder_pc, input_folder_res, filename)
% GENERATE_LABELS labels points belonging to different geometric primitives and writes the results to a file
%
% GENERATE_LABELS(input_folder_pc, input_folder_res, filename)
% [labels, types] = GENERATE_LABELS(input_folder_pc, input_folder_res, filename)
%
% Input:
% - input_folder_pc: path to folder containing the point cloud
% - input_folder_res: path to folder containing detection results
% - filename: filename of point cloud, without extension
%
% Output:
% - labels: list of labels (one per points)
% - types: list of primitive types (one per points)
% If no output is used, the results are written to [input_folder_res filename '_labels.dat'].

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

pc = pcread([input_folder_pc filename '.ply']);

labels = zeros(pc.Count, 1) - 1;
types = zeros(pc.Count, 1);
min_dist = .05 * ones(pc.Count, 1);

planes = load([input_folder_res filename '_planes.txt']);
for k=1:size(planes, 1)
    label = planes(k, 1);
    plane = planes(k, 2:end);
    d = dist_plane(pc, plane(1:3), plane(4));
    g = grad_plane(pc, plane(1:3), plane(4));
    idx = abs(d) < min_dist & abs(sum(pc.Normal .* g, 2)) > cosd(15);
    labels(idx) = label;
    types(idx) = 2; % plane
    min_dist(idx) = abs(d(idx));
end

spheres = load([input_folder_res filename '_spheres.txt']);
for k=1:size(spheres, 1)
    label = spheres(k, 1);
    sphere = spheres(k, 2:end);
    d = dist_sphere(pc, sphere(1:3), sphere(4));
    g = grad_sphere(pc, sphere(1:3), sphere(4));
    idx = abs(d) < min_dist & -sum(pc.Normal .* g, 2) > cosd(15);
    labels(idx) = label;
    types(idx) = 1; % sphere
    min_dist(idx) = abs(d(idx));
end
    
cyls = load([input_folder_res filename '_cylinders.txt']);
for k=1:size(cyls, 1)
    label = cyls(k, 1);
    cyl = cyls(k, 2:end);
    d = dist_cylinder(pc, cyl(1:3), cyl(4:6), cyl(7));
    g = grad_cylinder(pc, cyl(1:3), cyl(4:6), cyl(7));
    idx = abs(d) < min_dist & -sum(pc.Normal .* g, 2) > cosd(15);
    labels(idx) = label;
    types(idx) = 3; % cylinder
    min_dist(idx) = abs(d(idx));
end

cones = load([input_folder_res filename '_cones.txt']);
for k=1:size(cones, 1)
    label = cones(k, 1);
    cone = cones(k, 2:end);
    d = dist_cone(pc, cone(1:3), cone(4:6), cone(7));
    g = grad_cone(pc, cone(1:3), cone(4:6), cone(7));
    idx = abs(d) < min_dist & -sum(pc.Normal .* g, 2) > cosd(15);
    labels(idx) = label;
    types(idx) = 4; % cone
    min_dist(idx) = abs(d(idx));
end

if nargout < 1
    writematrix([types labels], [input_folder_res filename '_labels.txt'], 'Delimiter', ' ');
end