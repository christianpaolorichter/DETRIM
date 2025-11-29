function [x,y,s] = ...
    sim_2D_switching_brownian_trajectory(T,tau,D,locPrec)
% SIM_2D_SWITCHING_BROWNIAN_TRAJECTORY simulates a single-molecule
% trajectory undergoing 2-state switching between two different diffusion
% coefficients (D) and adding localization noise (locPrec).
%
% INPUTS:
%   T:       Total duration of the trajectory (number of frames).
%   tau:     [tau_0, tau_1] Mean lifetimes (in frames) of state 0 and state 1.
%   D:       [D_0, D_1] Diffusion coefficients (in px^2/frame) for state 0 and 1.
%   locPrec: Localization precision (standard deviation, in px) added as noise.
%
% OUTPUTS:
%   x:       Vector of observed X positions (with noise, in px).
%   y:       Vector of observed Y positions (with noise, in px).
%   s:       Vector of true kinetic states (0 or 1) for each frame.
%
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck

% --- 1. Initialization and Equilibrium Setup ---
pOcc = tau./sum(tau); % Calculate equilibrium state occupancy (probabilities of being in state 0 and state 1).
[x,y,s] = deal(nan(T,1)); % Pre-allocate memory for position (x, y) and true state (s) vectors.

x(1) = 0; y(1) = 0; % Initialize starting position at (0, 0).
if rand < pOcc(1) % Determine the initial state based on equilibrium occupancy probabilities.
    s(1) = 0; % Start in state 0 (e.g., Mobile)
else
    s(1) = 1; % Start in state 1 (e.g., Immobile)
end %if

% --- 2. Time-stepping Loop (Kinetic Switching and Diffusion) ---
for t = 2:T
    % --- Kinetic Switching Step ---
    % Calculate the probability of switching: 1 - exp(-dt / tau_i), where dt=1 frame.
    % If the random number is less than the switching probability, a switch occurs.
    if rand < 1-exp(-1/tau(s(t-1)+1))
        % Switch state: use mod(current state + 1, 2) to flip between 0 and 1.
        s(t) = mod(s(t-1)+1,2);
    else
        % No switch occurs: maintain the previous state.
        s(t) = s(t-1);
    end %if
    
    % --- Diffusion Step (Brownian Motion) ---
    % Draw the step size based on the current state's diffusion coefficient D.
    % Brownian step = randn * sqrt(2 * D * dt), where dt=1 frame.
    x(t) = x(t-1)+randn*sqrt(2*D(s(t)+1)); %[px] (Current position = previous position + random step)
    y(t) = y(t-1)+randn*sqrt(2*D(s(t)+1)); %[px]
end
    
% --- 3. Add Localization Precision Noise ---
% Add Gaussian noise with standard deviation 'locPrec' to all final positions.
x = x+randn(T,1)*locPrec; %[px] (Observed x position)
y = y+randn(T,1)*locPrec; %[px] (Observed y position)
end %fun