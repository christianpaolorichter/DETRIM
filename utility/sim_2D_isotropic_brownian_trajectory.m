function [x,y,t,dx,dy] = ...
    sim_2D_isotropic_brownian_trajectory(N,T,D,locPrec)
% SIM_2D_ISOTROPIC_BROWNIAN_TRAJECTORY Simulates N independent, two-dimensional (2D)  
%   isotropic Brownian trajectories, each of length T. It achieves this by calling a 1D
%   Brownian motion simulation function twice: once for the x-component and
%   once for the independent y-component.
%
%   Inputs:
%       N (integer): Number of independent trajectories to simulate.
%       T (integer): Length (number of time steps/frames) of each trajectory.
%       D (numeric): The 2D Diffusion coefficient [px^2/frame].
%       locPrec (numeric): Localization precision (measurement noise) [px].
%
%   Outputs:
%       x (matrix T x N): The simulated position data along the x-axis.
%       y (matrix T x N): The simulated position data along the y-axis.
%       t (matrix T x N): The time vector (frame number) corresponding to x/y.
%       dx (matrix T x N): The random step size (jump) at each time point in x.
%       dy (matrix T x N): The random step size (jump) at each time point in y.
%
%   Copyright (c) 2016-2025 Christian Paolo Richter
%   University of Osnabrueck

% Simulate the x-component of the 2D trajectory.
% For isotropic Brownian motion, the motion in X is independent of Y.
% The inputs N, T, D, and locPrec are passed directly.
% Outputs are x-position, time vector (t), and x-step size (dx).
[x,t,dx] = sim_1D_isotropic_brownian_trajectory(N,T,D,locPrec); %[px | frame]

% Simulate the y-component of the 2D trajectory.
% This is an entirely independent 1D simulation using the same parameters.
% Outputs are y-position, and y-step size (dy).
% The time vector (t) is calculated again by the 1D function but is discarded (~)
% since 't' from the first call is already correct and used.
[y,~,dy] = sim_1D_isotropic_brownian_trajectory(N,T,D,locPrec); %[px | frame]

end %fun