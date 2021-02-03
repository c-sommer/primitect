function im = create_overlay(im, pts, K)
% CREATE_OVERLAY projects a set of 3D points onto an image and displays it
%
% CREATE_OVERLAY(im, pts, K)
%
% Input:
% - im: image onto which points shall be projected
% - pts: Nx3 array containing coordinates of N 3D points
% - K: camera intrinsic matrix
%
% Output:
% - im: image with points projected onto it

% Published under GPL (v3+) License as part of PrimiTect project
% https://www.github.com/c-sommer/primitect/
% Copyright (c) 2019, Erik Bylow, Christiane Sommer.

% Assume camera is [I 0]
px = K*pts';
px = px./repmat(px(3,:),3,1);
px(1:2,:) = px(1:2,:) + 1;
px([1 2], :) = px([2 1], :);
px = int64(round(px));
px2 = px;
px3 = px;
px2(3,:) = px2(3,:) + 1;
px3(3,:) = px3(3,:) + 1;
ind = sub2ind(size(im), px(1,:), px(2,:), px(3,:));
ind2 = sub2ind(size(im), px2(1,:), px2(2,:), px2(3,:));
ind3 = sub2ind(size(im), px3(1,:), px3(2,:), px3(3,:));
im(ind) = 256;
im(ind2) = 0;
im(ind3) = 0;
end