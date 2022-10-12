# Learned-Compressed-Rank-Select-TALG22
This repository contains the code to reproduce the experiments in the [TALG22 paper](https://dl.acm.org/doi/pdf/10.1145/3524060):

> Antonio Boffa, Paolo Ferragina, and Giorgio Vinciguerra. 2022. A Learned Approach to Design Compressed Rank/Select Data Structures. ACM Trans. Algorithms 18, 3, Article 24 (October 2022), 28 pages. DOI:https://doi.org/10.1145/3524060

In brief, the paper follows a recent line of research on the so-called learned data structures. It provides a “learned” scheme for implementing a rank/select dictionary over compressed space. In particular, it introduces a novel lossless compressed storage scheme for the input dictionary which turns this problem into the one of approximating a set of points in the Cartesian plane via segments. The storage of the dictionary can be defined by means of a compressed encoding of these segments and the “errors” they do in approximating the input integers. Proper algorithms and data structures are then added to this compressed storage scheme to support fast rank and select operations.

A preliminary version of this work appeared in ([repo](https://github.com/aboffa/Learned-Rank-Select-ALENEX21)):

> Antonio Boffa, Paolo Ferragina, and Giorgio Vinciguerra. A "learned" approach to quicken and compress rank/select dictionaries. In Proceedings of the Symposium on Algorithm Engineering and Experiments (ALENEX). SIAM, 2020. DOI:https://doi.org/10.1137/1.9781611976472.4

This repo includes several new experiments: study of high-order compression of the aforesaid
approximation errors; a new improved hybrid data structure that combines our `la_vector` with existing rank/select dictionaries; more comprehensive experimental
evaluation of the `la_vector` that includes other recently proposed rank/select dictionary implementations (for example: sdsl::s18_vector ([repo](https://github.com/mudetz/s18_vector), [article](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=9281244)) and sdsl::rle_vector ([repo](https://github.com/vgteam/sdsl-lite)).

## Build and run

Clone this repo using  the flag `--recursive`. To run the experiments you need CMake 3.8+, and a compiler with support for C++17, Boost, and OpenMP.
To compile the executables, issue the following commands:

    ./lib/move_adapted_code.sh
    cmake . -B build -DCMAKE_BUILD_TYPE=Release
    cd build && make

The latter commands generates the executable for the benchmark (`my_benchmark`). To get the already manipulated datasets you can download them [here](https://drive.google.com/drive/folders/1K78tr9maRMPBhjx0Uo_SklogPY9wC7S5?usp=sharing).

The usage of `my_benchmark` is well explained in [this file](https://github.com/aboffa/Learned-Compressed-Rank-Select-TALG22/blob/main/include/arguments_parser.hpp).


The experiments can be run with the following script, which will populate a `result` directory with csv files:

    bash run_all.sh

The experiments require at least 32 GB of RAM and may take quite some time to finish.

## Tests

Since this benchmark deals with very different data structures implementations that have slightly different operations
there is a bunch of tests that check the actual correctness of every implementation and every wrapper.
They are in the directory `tests`. To run them:

    ./build/tests/my_tests

## Test environment

The code was tested on the following machine:

| Component | Specs                                     |
|-----------|-------------------------------------------|
| CPU       | Intel(R) Xeon(R) CPU E5-2407 v2 @ 2.40GHz |
| RAM       | 41 GB                                     |
| L1 cache  | 32 KB (data) 32 KB (instructions)         |
| L2 cache  | 256 KB                                    |
| L3 cache  | 1 MB                                      |
| OS        | Ubuntu 16.04.6 LTS                        |
| Compiler  | gcc 9.2.1                                 |
| CMake     | version 3.13.2                            |

## License

This project is licensed under the terms of the GNU General Public License v3.0.

If you use this code for your research, please cite:

```
@article{Boffa:2022,
author = {Boffa, Antonio and Ferragina, Paolo and Vinciguerra, Giorgio},
title = {A Learned Approach to Design Compressed Rank/Select Data Structures},
year = {2022},
issue_date = {July 2022},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {18},
number = {3},
issn = {1549-6325},
url = {https://doi.org/10.1145/3524060},
doi = {10.1145/3524060},
abstract = {We address the problem of designing, implementing, and experimenting with compressed data structures that support rank and select queries over a dictionary of integers. We shine a new light on this classical problem by showing a connection between the input integers and the geometry of a set of points in a Cartesian plane suitably derived from them. We then build upon some results in computational geometry to introduce the first compressed rank/select dictionary based on the idea of “learning” the distribution of such points via proper linear approximations (LA). We therefore call this novel data structure the la_vector. We prove time and space complexities of the la_vector in several scenarios: in the worst case, in the case of input distributions with finite mean and variance, and taking into account the kth order entropy of some of its building blocks. We also discuss improved hybrid data structures, namely, ones that suitably orchestrate known compressed rank/select dictionaries with the la_vector. We corroborate our theoretical results with a large set of experiments over datasets originating from a variety of applications (Web search, DNA sequencing, information retrieval, and natural language processing) and show that our approach provides new interesting space-time tradeoffs with respect to many well-established compressed rank/select dictionary implementations. In particular, we show that our select is the fastest, and our rank is on the space-time Pareto frontier.},
journal = {ACM Trans. Algorithms},
month = {oct},
articleno = {24},
numpages = {28},
keywords = {piecewise linear approximations, algorithm engineering, rank/select dictionaries, high order entropy, Compressed data structures}
}
```
