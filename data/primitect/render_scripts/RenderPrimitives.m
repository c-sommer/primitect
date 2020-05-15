%% setup

num_spheres = 20;
num_cylinders = 5;
num_cones = 10;
num_planes = 5;

% grid of dimensions 401 * 401
[x, y] = meshgrid(-2:.01:2,-2:.01:2);
sz = size(x);
x = x(:);
y = y(:);

num_points = sz(1) * sz(2);

%% create and sample spheres

k = 1;

params = cell(num_spheres + num_cylinders + num_cones + num_planes, 2);

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
    % replace z everywhere where sphere is in foreground
    z_max(idx) = z(idx);
    % replace normals and index as well
    i_max(idx) = k;
    % store current center and radius
    cr(k, :) = [c r];
    params{k, 1} = 1;
    params{k, 2} = [c r];
    k = k + 1;
end

while k <= num_spheres+num_cylinders
    c = rand(1, 3) .* [4 4 2] - [2 2 1];
    a = randn(1,3);
    a = a / norm(a);
    r = .05 + .95*rand(1);
    if abs(a(3)) > .8
        continue;
    end
    % compute z for according cylinder
    z = .1*x;
    sqr = r^2 - (x - c(1)).^2 - (y - c(2)).^2;
    scalp = (x - c(1)) * a(1) + (y - c(2)) * a(2);
    sqr = scalp.^2 + (1 - a(3)^2) * sqr;
    z(sqr > 0) = c(3) + (scalp(sqr > 0) * a(3) + sqrt(sqr(sqr > 0))) ./ (1 - a(3)^2);
    idx = z > z_max & z < 2;
    % replace z everywhere where cylinder is in foreground
    z_max(idx) = z(idx);
    % replace normals and index as well
    i_max(idx) = k;
    % store current center, axis and radius
    params{k, 1} = 3;
    params{k, 2} = [c a r];
    k = k + 1;
end

while k <= num_spheres+num_cylinders+num_cones
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
    params{k, 1} = 4;
    params{k, 2} = [c a asin(sin_theta)];
    k = k + 1;
end

while k <= num_spheres+num_cylinders+num_cones+num_planes
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
    % replace z everywhere where plane is in foreground
    z_max(idx) = z(idx);
    % replace normals and index as well
    i_max(idx) = k;
    % store current triangle
    params{k, 1} = 2;
    params{k, 2} = abc;
    k = k + 1;
end

params_new = cell(0, 2);

num_objects = k-1;

inliers = zeros(size(params, 1), 1);
for k=1:num_objects
    inliers(k) = sum(i_max == k);
    if inliers(k) < .02 * num_points
        inliers(k) = 0;
        i_max(i_max == k) = 0;
    else
        params_new{end+1, 1} = params{k, 1};
        params_new{end, 2} = params{k, 2};
    end
end

params = params_new;

% adapt labels accordingly

i_cr = [];
i_car = [];
i_cat = [];
i_nd = [];

z_max = .1*x;
i_max = zeros(size(x));
types = 2 * ones(size(x));
normal = ones(size(x,1), 3) .* [-.1 0 1];

for k=1:size(params, 1)
    switch params{k, 1}
        case 1 % sphere
            c = params{k, 2}(1:3);
            r = params{k, 2}(4);
            % compute z for according sphere
            z = .1*x;
            sqr = r^2 - (x - c(1)).^2 - (y - c(2)).^2;
            z(sqr >= 0) = c(3) + sqrt(sqr(sqr >= 0));
            idx = z > z_max;
            % replace z everywhere where sphere is in foreground
            z_max(idx) = z(idx);
            % replace normals and index as well
            normal(idx, :) = [x(idx) y(idx) z(idx)] - c;
            i_cr = [i_cr; k+1 c r];
        case 2 % plane
            abc = params{k, 2};
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
            % replace z everywhere where plane is in foreground
            z_max(idx) = z(idx);
            % replace normals and index as well
            n = n_l / norm(n_l);
            n = n * sign(n(3));
            normal(idx, :) = repmat(n, sum(idx), 1);
            i_nd = [i_nd; k+1 n -dot(n, abc(1,:))];
        case 3 % cylinder
            c = params{k, 2}(1:3);
            a = params{k, 2}(4:6);
            r = params{k, 2}(7);
            % compute z for according cylinder
            z = .1*x;
            sqr = r^2 - (x - c(1)).^2 - (y - c(2)).^2;
            scalp = (x - c(1)) * a(1) + (y - c(2)) * a(2);
            sqr = scalp.^2 + (1 - a(3)^2) * sqr;
            z(sqr > 0) = c(3) + (scalp(sqr > 0) * a(3) + sqrt(sqr(sqr > 0))) ./ (1 - a(3)^2);
            idx = z > z_max & z < 2;
            % replace z everywhere where cylinder is in foreground
            z_max(idx) = z(idx);
            % replace normals and index as well
            normal(idx, :) = [x(idx) y(idx) z(idx)] - c;
            normal(idx, :) = normal(idx, :) - (x(idx) * a(1) + y(idx) * a(2) + z(idx) * a(3) - dot(c, a)) * a;
            i_car = [i_car; k+1 c a r];
        case 4 % cone
            c = params{k, 2}(1:3);
            a = params{k, 2}(4:6);
            theta = params{k, 2}(7);
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
            % replace z everywhere where cone is in foreground
            z_max(idx) = z(idx);
            % replace normals and index as well
            normal(idx, :) = ([x(idx) y(idx) z(idx)] - c - h(idx) * a);
            normal(idx, :) = cos_theta * normal(idx, :) ./ sqrt(sum(normal(idx,:).^2, 2)) - sin_theta * a;
            i_cat = [i_cat; k+1 c a theta];
    end
    i_max(idx) = k;
    types(idx) = params{k, 1};
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

%% write clean object data to file

folder = '../all/';
stamp = num2str(round(1e6*(now-floor(now)))); %sprintf('%f', now);
% write to ply file
pcwrite(pc, [folder, 'noise000/', stamp, '.ply'], 'Encoding', 'binary');
% write to xyz file
pcwrite_xyz(pc, [folder, 'noise000/', stamp, '.xyz']);
% write labels to file
writematrix([types i_max+1], [folder, stamp, '_labels.txt'], 'Delimiter', '\t')
% write depth to png file
imwrite(.25*reshape(z_max, sz) + .5, [folder, stamp, '.png'], 'BitDepth', 16)
% write normals to png file
imwrite(.5*reshape(normal, [sz 3]) + .5, [folder, stamp, '_normals.png'])
% write labels to png file
imwrite(reshape(uint8(i_max)+1, sz), [folder, stamp, '_labels.png'])

% write parameters to files
n = [-.1 0 1] / norm([-.1 0 1]);
writematrix([[1 n 0]; i_nd], [folder, stamp, '_planes.txt'], 'Delimiter', '\t');
writematrix(i_cr, [folder, stamp, '_spheres.txt'], 'Delimiter', '\t');
writematrix(i_car, [folder, stamp, '_cylinders.txt'], 'Delimiter', '\t');
writematrix(i_cat, [folder, stamp, '_cones.txt'], 'Delimiter', '\t');

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