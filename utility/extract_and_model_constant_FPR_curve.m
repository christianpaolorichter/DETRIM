function fitresult = extract_and_model_constant_FPR_curve(win,D,FPR,targetFPR)
% EXTRACT_AND_MODEL_CONSTANT_FPR_CURVE interpolates the FPR surface and
% fits a power law (D = a * win^b) to the curve where FPR = targetFPR.
%
%   written by Christian Paolo Richter

%% Interpolate surface
[xData, yData, zData] = prepareSurfaceData(win,D,FPR);
ft = 'biharmonicinterp';
opts = fitoptions('Method','BiharmonicInterpolant');
opts.ExtrapolationMethod = 'biharmonic';
opts.Normalize = 'on';
fitresult = fit([xData, yData],zData,ft,opts);

% Create dense grid for interpolation
[win_,D_] = meshgrid(linspace(min(win),max(win),100),...
    linspace(min(D),max(D),100));
FPR_ = reshape(feval(fitresult,[win_(:),D_(:)]),100,100);

%% Find points with given FPR (D required for a constant targetFPR)
targetD = nan(1, 100);
for i = 1:100
    % Find the Diffusion Coefficient (D) where FPR_ equals targetFPR for a given winT
    % interp1 interpolates the D values based on the FPR_ surface
    % Note: Assumes FPR_ values are monotonically related to D for fixed winT
    targetD(i) = interp1(FPR_(:,i),D_(:,1),targetFPR);
end %for

%% Fit constant FPR curve to a Power Law (D = a * win^b)
take = not(isnan(targetD));
ft = fittype('power1'); % y = a*x^b
opts = fitoptions('Method','NonlinearLeastSquares');
opts.Display = 'Off';
fitresult = fit(colvec(win_(1,take)),...
    colvec(targetD(take)),ft,opts);
end %fun