DEtection of TRansient IMmobilization (DETRIM)
==========

DETRIM is a clustering algorithm implemented in MATLAB, designed specifically for analyzing Single-Molecule Localization Microscopy (SMLM) data to identify transient immobilization events.

---

## Package Components

This package consists of the following interdependent MATLAB functions:

| File Name | Purpose |
| :--- | :--- |
| `DETRIM.m` | **Main Entry Point.** Executes the hierarchical, multi-window search and iterative clustering. |
| `DETRIM_fwd_rev_cluster.m` | Performs the core clustering for a single time window, including **forward and reverse passes and fusion** via maximum core point overlap to mitigate temporal edge effects. |
| `DETRIM_group.m` | Implements the specialized DBSCAN-like grouping logic, ensuring cluster propagation only occurs through **Core Points**. |
| `DETRIM_pot_link.m` | Pre-filters the data using a k-d tree to find all purely spatial neighbors within `searchRad`. |
| `DETRIM_hard_thresh.m` | Applies the hard temporal filter to neighbors and **resolves simultaneous localizations** (multiple points in one frame) by selecting only the closest point. |
| `DETRIM_optimal_search.m` | Determines the optimal search radius and time search window given the minimal expected diffusion coefficient of particles and localization precision |

---

## Installation and Setup

### Prerequisites

1.  **MATLAB** installation (tested with 2024a).
2.  MATLAB's **Statistics and Machine Learning Toolbox**.

### Cloning the Repository

To get started, clone the repository to your local machine using Git Bash:

```bash
git clone git@github.com:christianpaolorichter/DETRIM.git
```

or using direct download:

1.  **Download ZIP:** Navigate to the repository page on GitHub. Click the **green `<> Code`** button and select **"Download ZIP."**
2.  **Unzip:** Extract the contents of the downloaded ZIP file to your desired project location (e.g., your MATLAB projects folder).

Then, add the main directory and all subfolders to your MATLAB environment:

```Matlab
% Run this command in MATLAB after cloning:
addpath(genpath('/path/to/DETRIM/'));
```

## Getting Started

The primary function to execute is `DETRIM.m`. This function requires only five input variables (`x`, `y`, `t`, `searchRad`, and `searchTime`) to run.

A toydata set was simulated (a single-molecule trajectory undergoing 2-state switching between mobile and immobile periods) and saved as localization_data.csv. In order to analyze this data using DETRIM:

```matlab
% --- 1. Load Localization Data ---
data = csvread('path\to\localization_data.csv');
x = data(:,1); %[px]     % X coordinates in pixels
y = data(:,2); %[px]     % Y coordinates in pixels
t = data(:,3); %[frame]  % Time frames

% --- 2. Calculate Spatial Search Radius (searchRad) ---
% searchRad is calculated based on localization precision (locPrec)
% and the desired probability (P) of finding an immobile particle's next localization.

pxSize = 107;   % [nm/px] Pixel size for converting units
locPrec = 25;   % [nm] Localization precision (as determined from your data)
P = 0.99;       % Expected True Positive Rate (99%)

% Formula ensures P% probability of re-detection for an immobile particle:
searchRad = 2 * locPrec * sqrt(-log(1 - P)) / pxSize; %[px] 

% --- 3. Set Temporal Search Window (searchTime) ---
% The window is set based on simulation-derived optimal values to balance 
% true-positive rate (>99%) and false-positive rate (<1%).

searchTime = 7; %[frame] 

% --- 4. Run DETRIM ---
% Use the 'verbose' flag (optional) for terminal feedback.
[ID,pntType] = DETRIM(...
    x, y, t, ...
    searchRad, ...
    searchTime, ...
    'verbose', true);
```
You should obtain the following Output (individual immobilization events are colored):
<p align="center">
  <img src="https://github.com/christianpaolorichter/DETRIM/blob/main/evaluation/result_DETRIM_localization_data.png?raw=true" alt=""/>
</p>

## Performance Evaluation (Simulation Results)

Extensive simulations were performed to assess the robustness, accuracy, and quantitative capability of DETRIM.

### I. True Positive Rate (Immobile Particles) and Kinetic Parameter Recovery

To assess the True Positive Rate (TPR) of the DETRIM algorithm, truly immobile particles ($D=0$) have been simulated with varying localization precision ($\epsilon$) ranging from 5 nm to 50 nm. For each condition, 1000 frames of continuous observation were generated acquiered with a frame rate of 32 ms. The lateral search radius (searchRad) was calculated to ensure a 99% probability of re-detecting the particle (see derivation in the manual). The immobile search was executed with a single pass, testing time search windows (searchTime) sequentially ranging from 2 to 10 frames. Each simulation condition was repeated 100 times to establish robust statistical accuracy.

<details>
<summary>See Evaluation Code</summary>

```matlab
% simulates immobile particles across a range of localization
% precisions and time search windows (winT) to calculate the True
% Positive Rate (TPR) of DETRIM.

pxSize = 107; %[nm] Pixel size for conversion
D = 0; %[nm^2/ms] immobile particle
locPrec = 5:5:50; %[nm] - Range of simulated localization precisions
lagTime = 32; %[ms]
searchTime = 2:10; % Range of minimum immobilization times (minT) to test
numIter = 30;

for idxPrec = numel(locPrec):-1:1
    for idxWin = numel(searchTime):-1:1
        for iter = numIter:-1:1 % Repeat each condition 30 times for statistical robustness
            TPR_(iter) = ...
                sim_and_apply_DETRIM(...
                D/pxSize^2*lagTime,...
                locPrec(idxPrec)/pxSize,...
                searchTime(idxWin),...
                'T',1000);
        end %for
        TPR(idxPrec,idxWin) = mean(TPR_);
    end %for
end

% --- Visualization ---
hFig = figure('color','w','Units','pixels',...
    'Position',set_figure_position(1,0.4,'center'),'visible','on');
hAx(1) = axes('Parent',hFig,'units','normalized','outerposition',[0,0,1,1],...
    'color','w','fontsize',12,'box','on','Nextplot','add','layer','top',...
    'xtick',searchTime,'ytick',locPrec);

% Plot the True Positive Rate (TPR) as a heatmap
IMG_plot(...
    (nanmask(TPR,TPR>0)),...
    'xdata',searchTime,...
    'ydata',locPrec,...
    'ColorMap',flipud(hot(22)),...
    'ColorBar','True Positive Rate [x100%]',...
    'ColorBarFontSize',12,...
    'ColorBarPosition','northoutside',...
    'EdgeLinestyle','-',...
    'EdgeColor',[0.7,0.7,0.7],...
    'hAx',hAx(1));

axis tight
xlabel('Time Search Window [frame]')
ylabel('Loc. Prec. [nm]')
clim([0.95 1])

disp('Mean TPR across all conditions:')
disp(mean(TPR,'all'))
```

</details>

* True positive rates (TPR) consistently **exceeded 99%** independent of the applied localization precision for search windows above 2 frames, demonstrating excellent performance in identifying immobilization events.

To assess DETRIM's utility in kinetic analysis, a stochastic bait-prey pull-down assay was simulated using 100 fixed bait sites over 1000 frames acquiered with a frame rate of 32 ms. The simulation modeled continuous association and dissociation with a mean bound lifetime ($\tau$) of 500 ms and a localization precision ($\epsilon$) of 25 nm. The lateral search radius (searchRad) was determined to guarantee a 99% probability of re-detecting bound particles (see derivation in the manual). The immobile search was executed using a single pass with a fixed time search window (searchTime) of 7 frames to ensure a >99% True Positive Rate (TPR).

<details>
<summary>See Evaluation Code</summary>

```matlab
% simulates a bait-prey pull-down assay to test DETRIM's 
% ability to successfully recover the underlying dissociation rate (tau) of the
% simulated binding events.

imgSize = 512; %[px] typical size for EMCCD sensor
pxSize = 107; %[nm] Pixel size for conversion
lagTime = 32; %[ms] - Time between frames
locPrec = 25; %[nm]
pAss = 0.01; % association probability per frame of prey binding to free bait
tau = 500; %[ms] average time of prey binding to bait (target value)
% Calculate dissociation probability (pDiss) per frame: pDiss = 1 - e^(-lagTime / tau)
pDiss = 1-exp(-lagTime/tau);
T = 1000; % Total number of frames (simulation duration)

% --- 1. Distribute Static Bait Positions ---
numBait = 100;
xPosBait = rand(numBait,1)*imgSize; %[px]
yPosBait = rand(numBait,1)*imgSize; %[px]

% --- 2. Run Stochastic Binding Kinetics Simulation ---
isOccupied = false(T,numBait); % state of bait occupation (0 = free, 1 = occupied by prey)
isOccupied(1,:) = rand(1,numBait) < pAss; % initialize state at t=1
for t = 2:T
    % Iterate through free bait: check for association
    for idxFree = rowvec(find(not(isOccupied(t-1,:))))
        if rand < pAss
            isOccupied(t,idxFree) = true;
        else
            isOccupied(t,idxFree) = false;
        end %if
    end %for
    
    % Iterate through occupied bait: check for dissociation
    for idxOccupied = rowvec(find(isOccupied(t-1,:)))
        if rand < pDiss
            isOccupied(t,idxOccupied) = false; % dissociation event
        else
            isOccupied(t,idxOccupied) = true;
        end %if
    end %for
end %for

% --- 3. Generate Localization Data from Occupied States ---
[t,idxBait] = find(isOccupied);
% Retrieve bait positions and apply finite localization precision noise
x = xPosBait(idxBait)+randn(numel(idxBait),1)*locPrec/pxSize; %[px]
y = yPosBait(idxBait)+randn(numel(idxBait),1)*locPrec/pxSize; %[px]

%% --- 4. DETRIM Clustering ---
P = 0.99;
% Calculate search radius based on target detection probability (P)
searchRad = 2*locPrec*sqrt(-log(1-P))/pxSize; %[px]
searchTime = 7; % Fixed minimum time window

clusterID = ...
    DETRIM(...
    x,y,t,...
    searchRad,...
    searchTime);

%% --- 5. Kinetic Analysis and Fitting ---
tauObs = [];
for idxCluster = max(clusterID):-1:2
    take = clusterID == idxCluster;
    % Calculate the observed duration (range in frames) for each cluster
    tauObs(idxCluster,1) = range(t(take));
end %for
tauObs(1) = []; % remove noise cluster (ID=1)

% Perform exponential fit (expfit)
% DETRIM only detects events of length >= searchTime (left-truncation).
% We fit the distribution of corrected times: (tauObs - searchTime + 1) * lagTime [ms]
[muhat,muci] = expfit((tauObs-searchTime+1)*lagTime);

% --- 6. Visualization ---
bins = 1:max(tauObs);
N = accumarray(tauObs,1); % Histogram of observed cluster durations

fprintf('Target Mean Lifetime (tau): %.1f ms\n', tau);
fprintf('Recovered Mean Lifetime (muhat): %.1f ms (95%% CI: %.1f - %.1f ms)\n', ...
    muhat, muci(1), muci(2));

figure; hold on

% Plot histogram of observed durations
plot(bins*lagTime,N,'k.','MarkerSize',10) 

% Plot theoretical fitted exponential decay curve (using the recovered muhat)
% The height is scaled by the total number of clusters (sum(N)), lagTime, and muhat
plot(bins*lagTime,...
    sum(N)/muhat*lagTime*exp(-((bins-searchTime+1)*lagTime/muhat)),'r','LineWidth',1.5)

title(sprintf('Total # detected Immobilization Events = %d',sum(N)))
xlabel('Observed Cluster Duration [ms]');
ylabel('Frequency (N)');
box on
grid on
axis tight

legend('Observed Cluster Count', ['Fitted Lifetime (\tau = ', ...
    num2str(round(muhat)), ' ms)'], 'Location', 'northeast');
```

</details>

* The recovered distribution of cluster durations was accurately fitted by an exponential decay model, and the **recovered mean bound lifetime ($\hat{\mu}$) closely matched the simulated target lifetime ($\tau$)**, proving the method's effectiveness for kinetic analysis.

<p align="center">
  <img src="https://github.com/christianpaolorichter/DETRIM/blob/main/evaluation/result_evaluate_DETRIM_TPR_and_tau.png?raw=true" alt=""/>
</p>

### II. False Positive Rate (Mobile Particles)

To assess DETRIM's robustness against misclassifying mobile events as immobile, mobile particles with diffusion coefficients ($D$) ranging from 10 $\text{nm}^2/\text{ms}$ to 200 $\text{nm}^2/\text{ms}$ were simulated. To mimic typical experimental conditions, a localization precision ($\epsilon$) of 25 nm and a frame lag time of 32 ms were used, with each condition repeated 100 times over 1000 frames for statistical robustness. The lateral search radius (searchRad) was set to ensure a 99% True Positive Rate (see derivation in the manual). The immobile search was executed using time search windows (searchTime) ranging from 2 to 20 frames.

<details>
<summary>See Evaluation Code</summary>

```matlab
% simulates mobile particles across a range of diffusion
% coefficients (D) and time search windows (winT) to characterize
% the False Positive Rate (FPR) of DETRIM.

pxSize = 107; %[nm] Pixel size for conversion
D = 10:10:200; %[nm^2/ms] - Range of simulated diffusion coefficients
locPrec = 35; %[nm] - Example localization precision used in the simulation
lagTime = 32; %[ms] - Time between frames
searchTime = 2:20; % Range of time search windows
numIter = 100; % Repeat each condition N times for statistical robustness

FPR = zeros(numel(D), numel(searchTime));

for idxD = numel(D):-1:1
    for idxWin = numel(searchTime):-1:1
        for iter = numIter:-1:1 
            FPR_(iter) = ...
                sim_and_apply_DETRIM(...
                D(idxD)/pxSize^2*lagTime,...
                locPrec/pxSize,...
                searchTime(idxWin),...
                'T',1000);
        end %for
        FPR(idxD,idxWin) = mean(FPR_);
    end %for
end %for

targetFPR = 0.01; % Target FPR threshold (1%)

% Fit the surface and extract the curve corresponding to the target FPR
fitresult = extract_and_model_constant_FPR_curve(searchTime,D,FPR,targetFPR);

%% Visualization
hFig = figure('color','w','Units','pixels',...
    'Position',set_figure_position(0.75,0.5,'center'),'visible','on');
hAx(1) = axes('Parent',hFig,'units','normalized','outerposition',[0,0,1,1],...
    'color','w','fontsize',12,'box','on','Nextplot','add','layer','top',...
    'ytick',searchTime,'xtick',D);

% Plot the False Positive Rate (FPR) as a heatmap (log10 scaled)
IMG_plot(...
    transpose(log10(nanmask(FPR,FPR>0))),...
    'ydata',searchTime,...
    'xdata',D,...
    'ColorMap',flipud(hot(22)),...
    'ColorBar','False Positive Rate (log_{10}) [x100%]',...
    'ColorBarFontSize',12,...
    'ColorBarPosition','northoutside',...
    'EdgeLinestyle','-',...
    'EdgeColor',[0.7,0.7,0.7],...
    'hAx',hAx(1));

axis tight
ylabel('Time Search Window [frame]')
xlabel('Diff. Coeff. [nm^2/ms]')
clim([-3 0]) % Sets color axis range from 1e-3 to 1e0 (0.1% to 100%)

hold on

% Plot the fitted curve corresponding to the targetFPR (1%)
plot(feval(fitresult,[searchTime searchTime(end)+0.5]),...
    [searchTime searchTime(end)+0.5],'b','linewidth',2)
axis([5.5 205.5 1.5 20.5])
```

</details>

* The time search window to achieve a false positive rate (FPR) below 1% at a fixed frame rate (here $\Delta t$ = 32 ms) depends on the diffusion coefficient ($D$) and the localization precision ($\epsilon$).

<p align="center">
  <img src="https://github.com/christianpaolorichter/DETRIM/blob/main/evaluation/result_evaluate_DETRIM_FPR.png?raw=true" alt=""/>
</p>

To set the optimal search radius and time search window given the minimal expected diffusion coefficient of particles and localization precision:
```matlab
pxSize = 107; %[nm] Pixel size for conversion
lagTime = 32; %[ms] Time interval between sequential image frames
D = 100; %[nm^2/ms] Diffusion Coefficient
locPrec = 35; %[nm] Localization Precision

% --- Optimal Search Parameter Calculation ---
[searchRad,searchTime] = ...
    DETRIM_optimal_search(...
    D/pxSize^2*lagTime,...
    locPrec/pxSize);

fprintf('Optimal lateral search radius [nm] = %.0f\n',searchRad*pxSize)
fprintf('Optimal time search window [frames] = %.0f\n',searchTime)    
```

## Citation
If you use this software/repository for your research, please cite:

Richter, C.P. (2025). *DETRIM* (v1.0.0). [Software]. Available from: https://github.com/christianpaolorichter/DETRIM