function d = min_dist(pc, planes, spheres, cylinders, cones)
% MIN_DIST computes the distance to the closest geometric primitive for each point
%
% d = MIN_DIST(pc, planes, spheres, cylindes, cones)
%
% Input:
% - pc: object of type PointCloud
% - planes: n_planes x 4 array containing plane parameters
% - spheresn n_spheres x 4 array containing sphere parameters
% - cylinders: n_cyl x 7 array containing cylinder parameters
% - cones: n_cones x 7 array containing cone parameters
%
% Output:
% - d: vector of distances to geometrically closest primitive (per point)

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Yumin Sun, Christiane Sommer.

N   = pc.Count;
Np  = size(planes, 1);
Ns  = size(spheres, 1);
Ncy = size(cylinders, 1);
Nco = size(cones, 1);
dist_plas = zeros(N, Np);
dist_sphs = zeros(N, Ns);
dist_cyls = zeros(N, Ncy);
dist_cons = zeros(N, Nco);

for p = 1:Np
    dist_plas(:,p) = dist_plane(pc,planes(p,1:3),planes(p,end));
end

for s = 1:Ns
    dist_sphs(:,s) = dist_sphere(pc,spheres(s,1:3),spheres(s,end));
end

for cy = 1:Ncy
    dist_cyls(:,cy) = dist_cylinder(pc,cylinders(cy,1:3),cylinders(cy,4:6),cylinders(cy,end));
end

for co = 1:Nco
    dist_cons(:,co) = dist_cone(pc,cones(co,1:3),cones(co,4:6),cones(co,end));
end

min_dist_plane    = min(abs(dist_plas), [], 2);
min_dist_sphere   = min(abs(dist_sphs), [], 2);
min_dist_cylinder = min(abs(dist_cyls), [], 2);
min_dist_cone     = min(abs(dist_cons), [], 2);
% Nx4
dist = [min_dist_plane, min_dist_sphere, min_dist_cylinder, min_dist_cone];
d = min(dist, [], 2);