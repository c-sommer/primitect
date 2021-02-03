function object_overlay(in_folder_pc, in_folder_res, in_folder_im, file_prefix, type)
% OBJECT_OVERLAY overlays a projection of 3D object inliers onto a given RGB image
%
% OBJECT_OVERLAY(in_folder_pc, in_folder_res, in_folder_im, file_prefix, type)
%
% Input:
% - in_folder_pc: string with folder containing the input ply file
% - in_folder_res: string with folder containing the result files
% - in_folder_im: string with folder containing the input RGB file
% - fileprefix: filename of point cloud without '.ply' (same as filename of RGB image without '_rgb.jpg')
%
% No output.

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Erik Bylow, Christiane Sommer.

fx = 525;
fy = 525;
cx = 319.5;
cy = 239.5;

K = [fx, 0 , cx; 0 fy cy; 0 0 1];

im = imread([in_folder_im, file_prefix, '_rgb.jpg']);
if strcmp(type, 'sphere')
    inliers = get_inliers_spheres(in_folder_pc, in_folder_res, file_prefix);
else
    if strcmp(type, 'cylinder')
        inliers = get_inliers_cylinder(in_folder_pc, in_folder_res, file_prefix);
    else
        if strcmp(type, 'cone')
            inliers = get_inliers_cone(in_folder_pc, in_folder_res, file_prefix);
        else
            error('Your specified type of primitive cannot be visualized')
        end
    end
end
im = create_overlay(im, inliers, K);
imshow(im);