%% spheres

for k=1:25
    RenderSpheres
    pause(rand(1))
end

%% cylinders

for k=1:25
    RenderCylinders
    pause(rand(1))
end

%% cones

for k=1:25
    RenderCones
    pause(rand(1))
end

%% planes

for k=1:25
    RenderPlanes
    pause(rand(1))
end

%% mixed

for k=1:100
    RenderPrimitives
    pause(rand(1))
end