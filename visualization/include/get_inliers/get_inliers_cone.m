function xyz = get_inliers_cone(in_folder_pc, in_folder_res, fileprefix)
% GET_INLIERS_CONE computes the points within a point cloud that are part of a given cone
%
% xyz = GET_INLIERS_CONE(in_folder_pc, in_folder_res, fileprefix)
%
% Input:
% - in_folder_pc: string with folder containing the input ply file
% - in_folder_res: string with folder containing the result files
% - fileprefix: filename of point cloud without '.ply'
%
% Output:
% - xyz: Nx3 array containing 3D coordinates of inliers
%
% See also: GET_INLIERS_SPHERE, GET_INLIERS_CYLINDER.

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Erik Bylow.

plyfile = [in_folder_pc fileprefix '.ply'];
cone_file = [in_folder_res fileprefix '_cones.txt'];

pc = pcread(plyfile);
cones = load(cone_file); %sx5
cone = cones(1, 2:end);

edges = linspace(0, .1, 1001);

dist = dist_cone(pc, cone(1:3), cone(4:6), cone(7));
grad = grad_cone(pc, cone(1:3), cone(4:6), cone(7));
idx = abs(dist) < .05 & -sum(pc.Normal .* grad, 2) > cosd(20);
xyz = pc.Location(idx, :);
end