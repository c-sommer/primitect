function pcshow_dist(pc, d, threshold)
% PCSHOW_DIST displays a point cloud where coloring is according to a given scalar field
%
% PCSHOW_DIST(pc, d)
% PCSHOW_DIST(pc, d, threshold)
%
% Input:
% - pc: pointCloud object with normals
% - d: distance or other scalar field to be displayed
% - threshold (optional): max value which shall still be resolved
%
% No output.

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

if nargin<3
    threshold = 0.05;
end

pcshow(pc.Location, min(abs(d), threshold))