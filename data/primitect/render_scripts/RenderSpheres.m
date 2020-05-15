%% setup

num_spheres = 30;

% grid of dimensions 401 * 401
[x, y] = meshgrid(-2:.01:2,-2:.01:2);
sz = size(x);
x = x(:);
y = y(:);

num_points = sz(1) * sz(2);

%% create and sample spheres

k = 1;
cr = zeros(num_spheres, 4);

z_max = .1*x;
i_max = zeros(size(x));

while k <= num_spheres
    c = rand(1, 3) .* [4 4 2] - [2 2 1];
    r = .05 + .95*rand(1);
    if c(3)+r < 0
        continue;
    end
    % compute z for according sphere
    z = .1*x;
    sqr = r^2 - (x - c(1)).^2 - (y - c(2)).^2;
    z(sqr >= 0) = c(3) + sqrt(sqr(sqr >= 0));
    idx = z > z_max;
    if all(~idx)
        continue;
    end
    % replace z everywhere where sphere is in foreground
    z_max(idx) = z(idx);
    % replace normals and index as well
    i_max(idx) = k;
    % store current center and radius
    cr(k, :) = [c r];
    k = k + 1;
end

cr_new = [];

inliers = zeros(size(cr, 1), 1);
for k=1:size(cr, 1)
    inliers(k) = sum(i_max == k);
    if inliers(k) < .02 * num_points
        inliers(k) = 0;
        i_max(i_max == k) = 0;
    else
        cr_new = [cr_new; cr(k, :)];
    end
end

cr = cr_new;

% adapt labels accordingly

z_max = .1*x;
i_max = zeros(size(x));
normal = ones(size(x,1), 3) .* [-.1 0 1];

for k=1:size(cr, 1)
    c = cr(k, 1:3);
    r = cr(k, 4);
    % compute z for according sphere
    z = .1*x;
    sqr = r^2 - (x - c(1)).^2 - (y - c(2)).^2;
    z(sqr >= 0) = c(3) + sqrt(sqr(sqr >= 0));
    idx = z > z_max;
    % replace z everywhere where sphere is in foreground
    z_max(idx) = z(idx);
    % replace normals and index as well
    normal(idx, :) = [x(idx) y(idx) z(idx)] - c;
    i_max(idx) = k;
end
normal = normal ./ sqrt(sum(normal.^2, 2));

pc = pointCloud(single([x y z_max]));
pc.Normal = single(normal);

%% visualization

subplot(121)
imagesc(reshape(z_max, sz))
axis image

subplot(122)
image(reshape(.5*(normal + 1), [sz 3]))
axis image

drawnow

%% write clean sphere data to file

folder = '../spheres/';
stamp = num2str(round(1e6*(now-floor(now)))); %sprintf('%f', now);
% write to ply file
pcwrite(pc, [folder, 'noise000/', stamp, '.ply'], 'Encoding', 'binary');
% write to xyz file
pcwrite_xyz(pc, [folder, 'noise000/', stamp, '.xyz']);
% write labels to file
types = ones(pc.Count, 1);
types(i_max==0) = 2;
writematrix([types i_max+1], [folder, stamp, '_labels.txt'], 'Delimiter', '\t')
% write depth to png file
imwrite(.25*reshape(z_max, sz) + .5, [folder, stamp, '.png'], 'BitDepth', 16)
% write normals to png file
imwrite(.5*reshape(normal, [sz 3]) + .5, [folder, stamp, '_normals.png'])
% write labels to png file
imwrite(reshape(uint8(i_max)+1, sz), [folder, stamp, '_labels.png'])

% write parameters to files
n = [-.1 0 1] / norm([-.1 0 1]);
writematrix([1 n 0], [folder, stamp, '_planes.txt'], 'Delimiter', '\t');
labels = (1:size(cr, 1)) + 1;
writematrix([labels', cr], [folder, stamp, '_spheres.txt'], 'Delimiter', '\t');
writematrix([], [folder, stamp, '_cylinders.txt']);
writematrix([], [folder, stamp, '_cones.txt']);

dlmwrite([folder 'files.txt'], stamp, 'delimiter', '', '-append')

%% add noise and write to files

noise_levels = .005:.005:.025;

for l=1:length(noise_levels)
    noise = noise_levels(l);
    pcn = pointCloud(pc.Location + noise_levels(l) * randn(pc.Count, 3));
    normals = pcnormals(pcn, 41);
    pcn.Normal = normals .* sign(sum(normals.*pc.Normal, 2));
    % write to ply file
    pcwrite(pcn, [folder, 'noise', sprintf('%03d', 1e3*noise), '/', stamp, '.ply'], 'Encoding', 'binary');
    % write to xyz file
    pcwrite_xyz(pcn, [folder, 'noise', sprintf('%03d', 1e3*noise), '/', stamp, '.xyz']);
end

%%

clear, clc