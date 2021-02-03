function grad = grad_cone(points, c, a, theta)
% GRAD_CONE
%
% grad = GRAD_CONE(points, c, a, theta)

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Christiane Sommer.

if isa(points, 'pointCloud')
    points = points.Location;
end

h = (points - c(:)') * a(:);
d_sq = sum((points - c(:)').^2, 2);
r = sqrt(d_sq - h.^2);
grad = zeros(size(points, 1), 3);
idx = h * cos(theta) + r * sin(theta) < 0;
grad(idx, :) = (c(:)' - points(idx, :)) ./ sqrt(d_sq(idx));
grad(~idx, :) = (sin(theta) + cos(theta) * h(~idx) ./ r(~idx)) * a - cos(theta) * (points(~idx, :) - c(:)') ./ r(~idx);