# ABM_SPT
Anchored Brownian motion for the analysis of single particle tracking data. Supporting code for the paper (Link to be added)

This source code does the following: 

1. Generates synthetic ABM data and initializes a MCMC chain.
2. Runs a tailor-made Gibbs sampler to generate MCMC samples from the posterior distribution.
3. Tracks MCMC samples and provides fit visualization.
4. Saves MCMC chain data for further exploration.

# Instructions (for use with user supplied data)
## Provide Data
The file `example_data.txt` contains data properly formatted for using the code. 
The format is as follows:
* The first line contains the time units and then the distance units, seperated by a comma.
* The second line encodes the order of the following data which is `m,t,wx,wy,hxx,hyy,hxy`.
* The remaining lines contain the data, seperated by commas in the order prescribed by the second line.

## Run code
If your data is in the form of `example_data.txt` then using the function `read_table(filename)` will read the data from this format, call the function `get_D_estimate` which initializes and runs the Gibbs sampler. `get_D_estimate` provides a posterior mean estimate of the diffusion coefficient, posterior variance estimate of the diffusion coefficient, and returns all the data from the MCMC run to a struct `chain`. The code that runs the Gibbs sampler is in the folder `2D_MP` which stands for "2-dimensional M trajectories." This path is added automatically by `get_D_estimate`.

# Instructions (for use with generated data)
## Run code
Within the folder `2D_MP`, the function `run_chainer(numTraj, len_t_global, flag_anchor)` generates syntethic data with `numTraj` particle trajectories which will have maximum length `len_t_global`. When `flag_anchor` is *TRUE* the analysis is conducted using anchored Brownian motion or and when *FALSE* the analysis is conducted using regular Brownian motion. `flag_anchor` is by default *TRUE*. `run_chainer` returns a struct `chain` containing all the data from the MCMC run. Each call of `run_chainer` generates new synthetic data using the file `generate_synthetic_data.`


