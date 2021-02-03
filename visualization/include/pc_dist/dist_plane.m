function dist = dist_plane(points, n, d)
% DIST_PLANE
%
% dist = DIST_PLANE(points, n, d)

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

if isa(points, 'pointCloud')
    points = points.Location;
end

dist = points * n(:) + d;