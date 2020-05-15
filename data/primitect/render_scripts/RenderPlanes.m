%% setup

num_planes = 20;

% grid of dimensions 401 * 401
[x, y] = meshgrid(-2:.01:2,-2:.01:2);
sz = size(x);
x = x(:);
y = y(:);

num_points = sz(1) * sz(2);

%% create and sample planes

k = 1;
nd = zeros(num_planes, 4);
abc_cell = cell(num_planes, 1);

z_max = .1*x;
i_max = zeros(size(x));

while k <= num_planes
    abc = rand(3,3) .* [4 4 2] - [2 2 0];
    b_a = abc(2,:) - abc(1,:);
    c_a = abc(3,:) - abc(1,:);
    n_l = cross(b_a, c_a);
    % compute z for according plane
    z = abc(1,3) - (n_l(1) * (x - abc(1,1)) + n_l(2) * (y - abc(1,2)))/n_l(3);
    p_a = [x y z] - abc(1,:);
    p_a = p_a ./ sqrt(sum(p_a.^2, 2));
    p_b = [x y z] - abc(2,:);
    p_b = p_b ./ sqrt(sum(p_b.^2, 2));
    p_c = [x y z] - abc(3,:);
    p_c = p_c ./ sqrt(sum(p_c.^2, 2));
    cos1 = sum(p_a .* p_b, 2);
    cos2 = sum(p_a .* p_c, 2);
    cos3 = sum(p_c .* p_b, 2);
    idx = z > z_max & (acos(cos1) + acos(cos2) + acos(cos3)) >= 2*pi*(1-eps);
    if all(~idx)
        continue;
    end
    % replace z everywhere where sphere is in foreground
    z_max(idx) = z(idx);
    % replace normals and index as well
    i_max(idx) = k;
    % store current triangle
    abc_cell{k} = abc;
    k = k + 1;
end

abc_cell_new = cell(0);

inliers = zeros(size(nd, 1), 1);
for k=1:size(nd, 1)
    inliers(k) = sum(i_max == k);
    if inliers(k) < .02 * num_points
        inliers(k) = 0;
        i_max(i_max == k) = 0;
    else
        abc_cell_new = [abc_cell_new; abc_cell{k}];
    end
end

abc_cell = abc_cell_new;

% adapt labels accordingly

z_max = .1*x;
i_max = zeros(size(x));
normal = ones(size(x,1), 3) .* [-.1 0 1];
nd = zeros(size(abc_cell, 1), 4);

for k=1:size(nd, 1)
    abc = abc_cell{k};
    b_a = abc(2,:) - abc(1,:);
    c_a = abc(3,:) - abc(1,:);
    n_l = cross(b_a, c_a);
    % compute z for according plane
    z = abc(1,3) - (n_l(1) * (x - abc(1,1)) + n_l(2) * (y - abc(1,2)))/n_l(3);
    p_a = [x y z] - abc(1,:);
    p_a = p_a ./ sqrt(sum(p_a.^2, 2));
    p_b = [x y z] - abc(2,:);
    p_b = p_b ./ sqrt(sum(p_b.^2, 2));
    p_c = [x y z] - abc(3,:);
    p_c = p_c ./ sqrt(sum(p_c.^2, 2));
    cos1 = sum(p_a .* p_b, 2);
    cos2 = sum(p_a .* p_c, 2);
    cos3 = sum(p_c .* p_b, 2);
    idx = z > z_max & (acos(cos1) + acos(cos2) + acos(cos3)) >= 2*pi*(1-eps);
    % replace z everywhere where sphere is in foreground
    z_max(idx) = z(idx);
    % replace normals and index as well
    n = n_l / norm(n_l);
    n = n * sign(n(3));
    nd(k, :) = [n -dot(n, abc(1,:))];
    normal(idx, :) = repmat(n, sum(idx), 1);
    i_max(idx) = k;
end

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

%% write clean plane data to file

folder = '../planes/';
stamp = num2str(round(1e6*(now-floor(now)))); %sprintf('%f', now);
% write to ply file
pcwrite(pc, [folder, 'noise000/', stamp, '.ply'], 'Encoding', 'binary');
% write to xyz file
pcwrite_xyz(pc, [folder, 'noise000/', stamp, '.xyz']);
% write labels to file
types = 2*ones(pc.Count, 1);
writematrix([types i_max+1], [folder, stamp, '_labels.txt'], 'Delimiter', '\t')
% write depth to png file
imwrite(.25*reshape(z_max, sz) + .5, [folder, stamp, '.png'], 'BitDepth', 16)
% write normals to png file
imwrite(.5*reshape(normal, [sz 3]) + .5, [folder, stamp, '_normals.png'])
% write labels to png file
imwrite(reshape(uint8(i_max)+1, sz), [folder, stamp, '_labels.png'])

% write parameters to files
n = [-.1 0 1] / norm([-.1 0 1]);
labels = (1:size(nd, 1)) + 1;
writematrix([[1 n 0]; [labels', nd]], [folder, stamp, '_planes.txt'], 'Delimiter', '\t');
writematrix([], [folder, stamp, '_spheres.txt']);
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