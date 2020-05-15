function [xyz] = range2points(D, K, zmax)
% RANGE2POINTS computes a point cloud from a given depth image
%
% xyz = RANGE2POINTS(D, K)
% xyz = RANGE2POINTS(D, K, zmax)
%
% Input:
% - D: depth image as returned by e.g. a Kinect
% - K: 3x3 intrinsic camera matrix
% - zmax: maximal depth that shall still be considered
%
% Output:
% - xyz: 3D points

% Christiane Sommer 05/2017

fx = K(1,1);
fy = K(2,2);
ux = K(1,3);
uy = K(2,3);

[n,m] = meshgrid(0:size(D,2)-1,0:size(D,1)-1);

if nargin<3
    idx = D > 0;
else
    idx = (D > 0) & (D < zmax);
end

x = D.*(n-ux)/fx;
y = D.*(m-uy)/fy;
z = D;

% xyz = [x(:),y(:),z(:)];
xyz = [x(idx),y(idx),z(idx)];