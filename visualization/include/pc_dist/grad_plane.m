function grad = grad_plane(points, n, d)
% GRAD_PLANE
%
% grad = GRAD_PLANE(points, n, d)

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

if isa(points, 'pointCloud')
    points = points.Location;
end

grad = repmat(n(:)', size(points,1), 1);