function dist = dist_cylinder(points, c, a, r)
% DIST_CYLINDER
%
% dist = DIST_CYLINDER(points, c, a, r)

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

if isa(points, 'pointCloud')
    points = points.Location;
end

dTa = (points - c(:)') * a(:);

dist = r - sqrt(sum((points - c(:)').^2, 2) - dTa.^2);