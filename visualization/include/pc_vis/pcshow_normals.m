function pcshow_normals(pc)
% PCSHOW_NORMALS displays a point cloud where coloring is chosen by normals
%
% PCSHOW_NORMALS(pc)
%
% Input:
% - pc: pointCloud object with normals
%
% No output.

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

pcshow(pc.Location, .5*(pc.Normal+1))