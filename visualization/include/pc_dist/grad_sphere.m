function grad = grad_sphere(points, c, r)
% GRAD_SPHERE
%
% grad = GRAD_SPHERE(points, c, r)

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

if isa(points, 'pointCloud')
    points = points.Location;
end

grad = (c(:)' - points) ./ sqrt(sum((points - c(:)').^2, 2));