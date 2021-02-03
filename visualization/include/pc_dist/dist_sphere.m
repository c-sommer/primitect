function dist = dist_sphere(points, c, r)
% DIST_SPHERE
%
% dist = DIST_SPHERE(points, c, r)

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

if isa(points, 'pointCloud')
    points = points.Location;
end

dist = r - sqrt(sum((points - c(:)').^2, 2));