# PrimiTect: Fast Continuous Hough Voting for Primitive Detection

This repository provides the code accompanying the paper
*PrimiTect: Fast Continuous Hough Voting for Primitive Detection*
by C. Sommer, Y. Sun, E. Bylow and D. Cremers,
presented at the International Conference on Robotics and Automation (ICRA) 2020.
A preprint can be found on [arXiv](https://arxiv.org/abs/2005.07457), and we will provide the DOI once available.

## Basic Usage

### Dependencies

The code depends on the following third-party libraries:

* Eigen (header-only)
* CLI (header-only)
* Ceres
* tinyply
* CGAL 4.14 (optional, only for comparing against Efficient RANSAC)

Eigen, Ceres and CLI are added to this repository as submodules, and from tinyply, we directly added the source files.

### Preparation

* Clone the repository to your computer including all submodules.
* Build Ceres in folder `thirdparty/build-ceres-solver/`.
* Compile the code using the `CMakeLists.txt` file:
    ```
    mkdir build
    cd build
    cmake ..
    make -j4
    cd ..
    ```

### Parameters

To display all input parameters for the respective executables, call `--help` in addition to the binary name.
In short, these are the most common parameters:

* `--folder`: the folder containing the input point cloud data
* `--in`: name of the input point clouds - must by in `.ply` format for our method, and in `.xyz` format for CGAL-based code
* `--results`: the folder in which results shall be written
* `--type`: primitive object type to be detected, the options are `plane`, `sphere`, `cylinder`, `cone`, `plane-sphere`, `plane-cyl`, `plane-cone` and `all`
* `--np` (default: 4096): number of input points to be used for voting stage, usually 2048-8192 is a good choice here
* `--settings` (default: all, only needed for our method): choose settings with which our code shall be run, `00` corresponds to the Drost* baseline implementation, while `39` corresponds to our proposed method, including vote spreading, cluster-by-averaging, and weighted averaging for candidate extraction

### Run Primitive Detection

#### Our proposed method

The main file, where our proposed method is run in its most basic form, is `cpp/detect_objects/main_voting_sample.cpp`.
The corresponding executable that is generated is `DetectObjects`.
To run it on a point cloud `test_data.ply` in folder `data`, `cd` into `cpp/detect_objects/`, and run
```
bin/DetectObjects --settings 39 --type all --in test_data.ply --folder ../../data --results .
```

#### Detection and refinement

The executables of type `SphereDetectRefine` etc. run a subsequent refinement after primitive detection.
Note that this is not part of our proposed method, but we include it to provide some way to get more accurate estimates.
It uses Ceres with a robust loss function.

#### Refinement of pre-detected primitives

The executables of type `OptSphere` etc. do only refinement, i.e. they require as input a list of primitive parameters, e.g. obtained from our detection method.
Note that these are also only provided for completeness and not part of the method we propose in the paper.

#### Primitive tracking

In the folder `cpp/track_object/`, you can find the code we used to track a single geometric primitive in a Redwood sequence, as shown in our supplementary video.
The folder contains code for detection with both our method (including the baseline) and the CGAL Efficient RANSAC implementation.
If you want the Efficient RANSAC files to build, you will have to uncomment the corresponding lines in `cpp/CMakeLists.txt` and `cpp/track_object/CMakeLists.txt`.

#### Efficient RANSAC for comparison

Finally, we include the code to run Efficient RANSAC on a point cloud in `.xyz` format, and convert results into our format, for easier comparison.
The code can be found in `cpp/cgal_test/` and is an adaptation from the CGAL examples on Efficient RANSAC.
If you want these files to build automatically, you will have to uncomment the corresponding lines in `cpp/CMakeLists.txt`.

## Data

### Redwood data (*primitive-rw*)

We provide the scripts we were using to obtain the data from the [Redwood dataset](http://redwood-data.org/3dscan/dataset.html) in the folder `data/redwood_lod/`.
You need to download the `.zip` folders containing the according sequences first, but then can run `extract_first_frames.sh`.

The Matlab script `RedwoodToPly.m` can then convert the depth images to `.ply` data, as required for our method.

### Synthetic data (*primitect*)

The Matlab script `RenderSyntheticAll.m` in `data/primitect/` will generate random point clouds of the same type as those we were using for our evaluation, and refer to as `primitect` dataset in the paper.


## License and Publication

Our code is released under the GPL (v3+) license, for more details please see the `LICENSE` file.
Also note the different licenses of the submodules in the folder `thirdparty`.

Please cite our paper when using the code in a scientific project. You can copy-paste the following BibTex entry:

```
@inproceedings{sommer2020,
    title   = {PrimiTect: Fast Continuous Hough Voting for Primitive Detection},
    author  = {Sommer, Christiane and Sun, Yumin and Bylow, Erik and Cremers, Daniel},
    booktitle = {IEEE International Conference on Robotics and Automation (ICRA)},
    year    = {2020}
}
```
