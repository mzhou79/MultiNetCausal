# Algorithm 1

Code updated to follow Algorithm 1 in the manuscript. Treatment assignment vectors are generated outside the Gibbs loop and held fixed within each chain.

With B = 1, estimates were similar to those from the manuscript's simulation implementation in our numerical checks, although the treatment-allocation update schemes differ.

examples/algorithm1.R illustrates the complete workflow using synthetic data: parameter estimation, effect estimation, and parametric bootstrap.