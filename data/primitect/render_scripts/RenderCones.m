%% setup

num_cones = 30;

% grid of dimensions 401 * 401
[x, y] = meshgrid(-2:.01:2,-2:.01:2);
sz = size(x);
x = x(:);
y = y(:);

num_points = sz(1) * sz(2);

%% create and sample cones

k = 1;
cat = zeros(num_cones, 7);

z_max = .1*x;
i_max = zeros(size(x));

while k <= num_cones
    c = rand(1, 3) .* [4 4 2] - [2 2 0];
    a = randn(1,3) .* [1 1 2];
    a = a / norm(a);
    sin_theta = sind(10) + (sind(45) - sind(10)) * rand(1);
    cos_theta = sqrt(1 - sin_theta^2);
    if cos_theta < a(3)
        continue;
    end
    % compute z for according cone
    z = .1*x;
    sqr = (x - c(1)).^2 + (y - c(2)).^2;
    scalp = (x - c(1)) * a(1) + (y - c(2)) * a(2);
    sqr = scalp.^2 * a(3)^2 - (a(3)^2 - cos_theta^2) * (scalp.^2 - cos_theta^2 * sqr);
    z(sqr > 0) = c(3) + (-scalp(sqr > 0) * a(3) - sqrt(sqr(sqr > 0))) ./ (a(3)^2 - cos_theta^2);
    h = ones(size(z));
    h(sqr > 0) = scalp(sqr > 0) + (z(sqr > 0) - c(3)) * a(3);
    idx = z > z_max & z < 2 & h > 0 & h < .8*sqrt(2*cos_theta^2/sin_theta);
    if all(~idx)
        continue;
    end
    % replace z everywhere where cone is in foreground
    z_max(idx) = z(idx);
    % replace normals and index as well
    i_max(idx) = k;
    % store current center, axis and angle
    cat(k, :) = [c a asin(sin_theta)];
    k = k + 1;
end

cat_new = [];

inliers = zeros(size(cat, 1), 1);
for k=1:size(cat, 1)
    inliers(k) = sum(i_max == k);
    if inliers(k) < .02 * num_points
        inliers(k) = 0;
        i_max(i_max == k) = 0;
    else
        cat_new = [cat_new; cat(k, :)];
    end
end

cat = cat_new;

% adapt labels accordingly

z_max = .1*x;
i_max = zeros(size(x));
normal = ones(size(x,1), 3) .* [-.1 0 1];

for k=1:size(cat, 1)
    c = cat(k, 1:3);
    a = cat(k, 4:6);
    theta = cat(k, 7);
    cos_theta = cos(theta);
    % compute z for according cone
    z = .1*x;
    sqr = (x - c(1)).^2 + (y - c(2)).^2;
    scalp = (x - c(1)) * a(1) + (y - c(2)) * a(2);
    sqr = scalp.^2 * a(3)^2 - (a(3)^2 - cos_theta^2) * (scalp.^2 - cos_theta^2 * sqr);
    z(sqr > 0) = c(3) + (-scalp(sqr > 0) * a(3) - sqrt(sqr(sqr > 0))) ./ (a(3)^2 - cos_theta^2);
    h = ones(size(z));
    h(sqr > 0) = scalp(sqr > 0) + (z(sqr > 0) - c(3)) * a(3);
    idx = z > z_max & z < 2 & h > 0 & h < .8*sqrt(2*cos_theta^2/sin_theta);
    % replace z everywhere where cylinder is in foreground
    z_max(idx) = z(idx);
    % replace normals and index as well
    normal(idx, :) = ([x(idx) y(idx) z(idx)] - c - h(idx) * a);
    normal(idx, :) = cos_theta * normal(idx, :) ./ sqrt(sum(normal(idx,:).^2, 2)) - sin_theta * a;
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

%% write clean cone data to file

folder = '../cones/';
stamp = num2str(round(1e6*(now-floor(now)))); %sprintf('%f', now);
% write to ply file
pcwrite(pc, [folder, 'noise000/', stamp, '.ply'], 'Encoding', 'binary');
% write to xyz file
pcwrite_xyz(pc, [folder, 'noise000/', stamp, '.xyz']);
% write labels to file
types = 4*ones(pc.Count, 1);
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
writematrix([], [folder, stamp, '_spheres.txt']);
writematrix([], [folder, stamp, '_cylinders.txt']);
labels = (1:size(cat, 1)) + 1;
writematrix([labels', cat], [folder, stamp, '_cones.txt'], 'Delimiter', '\t');

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