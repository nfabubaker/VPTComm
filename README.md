<img src="https://drive.google.com/uc?export=view&id=1nF7gGV1YjqOp6mV7twIcH0_ES-KLd3TX" alt="logo" width="600" height="200">

## Overview

**VPTComm** is a library that optimizes MPI point-to-point communication by reducing message counts using a k-ary n-cube virtual process topology and store-and-forward techniques. This library helps mitigate high latency in MPI applications by significantly reducing the number of messages exchanged between processes.

Written in plain C and dependent only on MPI, **VPTComm** provides a simple and efficient way to replace traditional point-to-point communication with a more optimized, low-latency alternative.

## Features

- **Optimized Communication**: Reduces message counts using a virtual k-ary n-cube topology.
- **Store-and-Forward**: Minimizes communication by employing store-and-forward techniques.
- **Compatibility**: Fully compatible with MPI, no additional dependencies are required.
- **Low Latency**: Achieves logarithmic message reduction for improved performance.
- **Tested**: Includes built-in tests to verify functionality and ensure reliable communication.

## Installation

To install **VPTComm**, clone this repository and follow the steps below:

### Prerequisites

- MPI (e.g., OpenMPI, MPICH)
- A C compiler (e.g., GCC)

### Installation Steps

1. Clone the repository:
    ```bash
    git clone https://github.com/nfabubaker/vptcomm.git
    cd vptcomm
    ```

2. Build the library and tests:
    ```bash
    mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=/path/to/dir ..
    make && make run_test && make install
    ```

## How to Use

1. **Include the library in your MPI application**:
    Add the following line to include **VPTComm**:
    ```c
    #include "vptcomm.h"
    ```
**Important note:** the current data type is set to "double", you should change it according to your application. Modify lines 13 & 16 in vptcomm.h

The library can be used as follows
```c

STFW_init(int _nInstances, unsigned char dual_enabled);
/* This call initializes the VPTComm Library. It is usually called just after MPI_Init in your application.
The _nInstances parameter tells the library how many P2P communication instances you would like to call in your application.
The dual_enabled paramteter tells the library if there are dual P2P communication calls (such as in an Expand/Reduce scenario)*/
.
.
.

STFW_init_instance(int _instance_id, int vpt_ndims, int *vpt_dsizes, int npsend, int *sendlist, int nprecv, int *recvlist, int *ssend, vptcomm_real_t **sendp, int *srecv, vptcomm_real_t **recvp, int *topomap);
/*This is the call that replaces your MPI_Send/Recv, MPI_Alltoall_v, or MPI_Neighbor_alltoallv calls.
It has similar parameters: */
.
.
.

STFW_Comm(int _instance_id);
/*The call starts the communication. If your application does not repeast the P2P communication more than once,
this has to be called directly after the the STFW_init_instance call. If your application repeats the P2P communication pattern (even if the data change), then the STFW_init_instance can be called once, and the STFW_Comm can be called inside the repetition loop.*/
.
.
.

STFW_finalize();
/*Finalizes the VPTComm usage, frees the allocated memory and distroys all relevant structs*/
```


2. **Link the library during compilation**:
    When building your application, link to the **vptcomm** library:
    ```bash
    mpicc -o my_application my_application.c -L. -lvptcomm
    ```

3. **Run your application**:
    Use `mpirun` to run the program with multiple processes:
    ```bash
    mpirun -np 8 ./my_application
    ```


## Developers

- [**Nabil Abubaker**--ETH Zurich](https://github.com/nfabubaker) and [**R. Oguz Selvitopi** -- LBNL](https://github.com/roguzsel)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
