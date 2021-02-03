function dist = dist_cone(points, c, a, theta)
% DIST_CONE
%
% dist = DIST_CONE(points, c, a, theta)

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

if isa(points, 'pointCloud')
    points = points.Location;
end

h = (points - c(:)') * a(:);
d_sq = sum((points - c(:)').^2, 2);
r = sqrt(d_sq - h.^2);
dist = zeros(size(points, 1), 1);
idx = h * cos(theta) + r * sin(theta) < 0;
dist(idx) = -sqrt(d_sq(idx));
dist(~idx) = h(~idx) * sin(theta) - r(~idx) * cos(theta);