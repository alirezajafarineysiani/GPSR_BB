# C and Python Implementation of Gradient Projection for Sparse Reconstruction (GPSR) Algorithm

## Overview

Welcome to the GPSR_BB repository, where you can find implementations of the "Gradient Projection for Sparse Reconstruction" (GPSR_BB) algorithm in Python and C. Originally developed in MATLAB by Mário Figueiredo, Robert D. Nowak, and Stephen J. Wright, this project aims to provide accessible versions of the algorithm in both programming languages.

## Sparse Recovery Algorithm

Sparse recovery addresses the fundamental problem in signal processing of reconstructing a sparse signal from limited measurements. The GPSR_BB algorithm is renowned for its effectiveness in efficiently solving this problem.

The problem is formulated as:

$$ \underset{\mathbf{x}}{\text{minimize}} \quad \lVert\mathbf{x}\rVert_0 \text{subject to} \quad \mathbf{A}\mathbf{x} = \mathbf{b}$$

Due to the NP-hard nature of the above problem, approaches based on the relaxation of the $\ell_0$-norm by the $\ell_1$-norm have been proposed.
GPSR solves the lasso problem, which is defined as follows:

$$ \underset{\mathbf{x}}{\text{minimize}} \quad \tau \lVert\mathbf{x}\rVert_1 + \lVert\mathbf{A}\mathbf{x} - \mathbf{b}\rVert_2^2$$

## Getting Started

To use these implementations, follow the steps below:

### Prerequisites

#### OpenBLAS Library (for C Implementation)

To build and run the C implementation, you need the OpenBLAS library. Here's how you can install it:

1. **Download OpenBLAS:**
   ```bash
   wget https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.27/OpenBLAS-0.3.27.tar.gz
   ```
   
   

3. **Extract the Archive:**
   ```bash
   tar -xzvf OpenBLAS-0.3.27.tar.gz
   cd OpenBLAS-0.3.27
   ```
   

4. **Configure OpenBLAS:**
   ```bash
   make TARGET=your_architecture
   ```
   
   Replace \`your_architecture\` with your specific architecture (e.g., HASWELL). Supported architectures can be found in the \`TARGET\` folder.

5. **Build and Install OpenBLAS:**
   ```bash
   make
   sudo make install
   ```
   This command installs the library system-wide.

### Cloning the Repository

   Clone the GPSR_BB repository using Git:
   ```bash
   git clone https://github.com/alirezajafarineysiani/GPSR_BB.git
   ```

### Building the C Implementation

To build the C implementation (example assuming gcc and OpenBLAS installed):
```bash
gcc main.c -lopenblas -L/opt/OpenBLAS/lib -I/opt/OpenBLAS/include -o gpsr_bb
```

Replace \`/opt/OpenBLAS/lib\` and \`/opt/OpenBLAS/include\` with the actual paths where OpenBLAS is installed.

### Using the Python Implementation

The Python implementation is provided as \`gpsr_bb.py\`. You can use it directly in your Python environment.

## Contributors

- [Alireza Jafari-Neysiani](https://github.com/alirezajafarineysiani)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
EOF
