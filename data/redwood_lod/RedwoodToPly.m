%% init

folder = 'cones/';
stamps = load([folder, 'stamps.dat']);

K = [525 0 319.5; 0 525 239.5; 0 0 1];

%% go through point clouds, compute normals and save

for k=1:length(stamps)
    stamp_k = sprintf('%05i', stamps(k));
    disp(stamp_k)
    D = imread([folder stamp_k '_depth.png']);
    xyz = range2points(double(D)/1000, K);
    pc = pointCloud(single(xyz));
    n = pcnormals(pc, 100); % approximately 10 pixel neighborhood
    pc.Normal = -n .* sign(sum(n.*xyz, 2));
    pcwrite(pc, [folder 'full/' stamp_k '.ply'], 'Encoding', 'binary');
    pcwrite_xyz(pc, [folder 'full/' stamp_k '.xyz']);
    pc = select(pc, find(pc.Location(:, 3)<2));
    pcwrite(pc, [folder 'cut2m/' stamp_k '.ply'], 'Encoding', 'binary');
    pcwrite_xyz(pc, [folder 'cut2m/' stamp_k '.xyz']);
end