function [x,t,dx] = ...
    sim_1D_isotropic_brownian_trajectory(N,T,D,locPrec)
% SIM_1D_ISOTROPIC_BROWNIAN_TRAJECTORY simulates N independent, one-dimensional (1D) Brownian
%   trajectories, each of length T. It calculates the random step size
%   (dx) based on the diffusion coefficient (D) and then integrates these
%   steps to find the position (x). It optionally adds localization noise.
%
%   Inputs:
%       N (integer): Number of independent trajectories to simulate.
%       T (integer): Length (number of time steps/frames) of each trajectory.
%       D (numeric): The Diffusion coefficient [px^2/frame].
%       locPrec (numeric): Localization precision (standard deviation of
%                          measurement noise) [px]. If 0, no noise is added.
%
%   Outputs:
%       x (matrix T x N): The simulated position (trajectory) data.
%       t (matrix T x N): The time vector (frame number) corresponding to x.
%       dx (matrix T x N): The random step size (jump) at each time point.
%
%   Copyright (c) 2014-2025 Christian Paolo Richter
%   University of Osnabrueck

% --- Input Validation ---
% Initialize the input parser object for validating inputs.
objParser = inputParser;

% Add validation rules for required inputs:
% N: Must be numeric, an integer (rem(x,1)==0), and >= 1.
objParser.addRequired('N',@(x) isnumeric(x) && rem(x,1) == 0 && x >= 1)
% T: Must be numeric, an integer (rem(x,1)==0), and >= 1.
objParser.addRequired('T',@(x) isnumeric(x) && rem(x,1) == 0 && x >= 1)
% D: Must be non-negative (Diffusion coefficient).
objParser.addRequired('D',@(x) x >= 0)
% locPrec: Must be numeric and non-negative (Localization precision).
objParser.addRequired('locPrec',@(x) isnumeric(x) && x >= 0)

% Execute the validation of the inputs.
objParser.parse(N,T,D,locPrec);

%% Simulation Core

% Calculate the random step size (dx) for T time steps and N trajectories.
% For a 1D Brownian motion: $\text{step} = \mathcal{N}(0, \sigma^2)$, where $\sigma^2 = 2D \Delta t$.
% Assuming $\Delta t = 1$ (frame duration): $\text{step} = \mathcal{N}(0, 2D)$.
% 'randn(T,N)' generates T x N matrix of random numbers from $\mathcal{N}(0, 1)$.
% Multiplying by 'sqrt(2*D)' scales the standard deviation to $\sqrt{2D}$.
dx = randn(T,N)*sqrt(2*D); %for 1dim jump process

% Calculate the position 'x' by integrating the step sizes 'dx'.
% 'cumsum(dx,1)' performs a cumulative sum down the first dimension (time, T).
% This integrates the random steps to get the trajectory position over time.
x = cumsum(dx,1); %[px]

% Check if localization noise (measurement error) should be added.
if locPrec > 0
    % Add localization noise: The observed position is the true position
    % plus noise drawn from $\mathcal{N}(0, \text{locPrec}^2)$.
    % 'randn(T,N)*locPrec' generates T x N noise matrix with $\sigma = \text{locPrec}$.
    x = x+randn(T,N)*locPrec;
end %if

% Create the time vector 't'.
% 'colvec(1:T)' (assumed function) creates a column vector of time steps [1; 2; ... T].
% Multiplying by 'ones(1,N)' replicates this column vector N times to match
% the dimensions of 'x' and 'dx' (T x N).
t = colvec(1:T)*ones(1,N);

end %fun