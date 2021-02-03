function xyz = get_inliers_sphere(in_folder_pc, in_folder_res, fileprefix)
% GET_INLIERS_SPHERE computes the points within a point cloud that are part of a given sphere
%
% xyz = GET_INLIERS_SPHERE(in_folder_pc, in_folder_res, fileprefix)
%
% Input:
% - in_folder_pc: string with folder containing the input ply file
% - in_folder_res: string with folder containing the result files
% - fileprefix: filename of point cloud without '.ply'
%
% Output:
% - xyz: Nx3 array containing 3D coordinates of inliers
%
% See also: GET_INLIERS_CYLINDER, GET_INLIERS_CONE.

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Erik Bylow.

plyfile = [in_folder_pc fileprefix '.ply'];
sph_file = [in_folder_res fileprefix '_spheres.txt'];

pc = pcread(plyfile);
sphs = load(sph_file); %sx5
sphere = sphs(1, 2:end);

edges = linspace(0, .1, 1001);

dist = dist_sphere(pc, sphere(1:3), sphere(4));
grad = grad_sphere(pc, sphere(1:3), sphere(4));
idx = abs(dist) < .05 & -sum(pc.Normal .* grad, 2) > cosd(20);
xyz = pc.Location(idx, :);
end