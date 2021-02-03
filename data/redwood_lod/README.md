We provide the scripts we were using to obtain the data from the [Redwood dataset](http://redwood-data.org/3dscan/dataset.html) here.
You need to download the `.zip` folders containing the according sequences first, and then can run `extract_first_frames.sh`.

The Matlab script `RedwoodToPly.m` converts the depth images to `.ply` data, as required for our method.

The files `cylinders/cut2m/01220.ply` and `cylinders/cut2m/01220.ply` are example files that can be used to test the algorithm without downloading all data from the redwood website..