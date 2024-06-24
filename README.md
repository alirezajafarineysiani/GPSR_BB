# C and Python implementation of Gradient Projection for Sparse Reconstruction (GPSR_BB)

## Overview

This repository contains implementations of the "Gradient Projection for Sparse Reconstruction" (GPSR_BB) algorithm in both Python and C. Originally developed in MATLAB by Mário Figueiredo, Robert D. Nowak, and Stephen J. Wright, this project aims to provide convenient and accessible versions of the algorithm in both programming languages.

## Sparse Recovery Algorithm

Sparse recovery is a fundamental problem in signal processing and related fields, where the goal is to accurately reconstruct a sparse signal from a limited set of measurements. The Gradient Projection for Sparse Reconstruction (GPSR_BB) algorithm is a popular method for achieving this task.

The problem is formulated as:

$$ \underset{\textbf{x}}{\text{minimize}} \quad \lVert \textbf{x} \rVert_0 \quad
\text{subject to} \quad \textbf{A}\textbf{x} = \textbf{b} $$

## Getting Started

To use these implementations, follow these steps:

### Cloning the Repository

1. Clone the repository:

   ```bash
   git clone https://github.com/alirezajafarineysiani/GPSR_BB.git
