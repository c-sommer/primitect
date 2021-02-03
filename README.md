# PrimiTect: Fast Continuous Hough Voting for Primitive Detection

This repository provides the code accompanying the paper
[*PrimiTect: Fast Continuous Hough Voting for Primitive Detection*](https://doi.org/10.1109/ICRA40945.2020.9196988)
by C. Sommer, Y. Sun, E. Bylow and D. Cremers,
presented at the International Conference on Robotics and Automation (ICRA) 2020.
A preprint can be found on [arXiv](https://arxiv.org/abs/2005.07457).

## Basic Usage

We provide detailed information on how to build, compile and run our code in the [cpp](cpp/) folder containing the code.

## Data

See the folder [data](data/) for information and scripts to generate the data we used in our paper.

## Visualization of Results

The figures in our paper can be reproduced using the Matlab visualization files provided in [visualization](visualization/).
We also provide example data (in `data/`) together with example results (in `ex_results`) which can be directly used to test the visualization functions.
Simply `cd` into `visualization/`, open Matlab and run the script `VisualizeExampleResults`.

## License and Publication

Our code is released under the GPL (v3+) license, for more details please see the `LICENSE` file.
Also note the different licenses of the submodules in the folder `thirdparty`.

Please cite our paper when using the code in a scientific project. You can copy-paste the following BibTex entry:

```
@inproceedings{sommer2020,
    title   = {PrimiTect: Fast Continuous Hough Voting for Primitive Detection},
    author  = {Sommer, Christiane and Sun, Yumin and Bylow, Erik and Cremers, Daniel},
    booktitle = {IEEE International Conference on Robotics and Automation (ICRA)},
    year    = {2020},
    doi = {10.1109/ICRA40945.2020.9196988}
}
```