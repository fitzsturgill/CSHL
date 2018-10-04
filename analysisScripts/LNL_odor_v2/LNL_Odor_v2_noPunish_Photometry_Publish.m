
function TE = LNL_Odor_v2_noPunish_Photometry_Publish(DB, animal, reloadPhotometry, 
% LNL_Odor_v2_noPunish_Photometry_Publish
% redo pupil, whisk, wheel, reject?


%% for deconvolution
tau = 2; % time constant of exponential decay kernel
kDuration = 3; % duration of kernel
Fs = 20;
kt = (1:kDuration*20)*(1/Fs) - (1/Fs); 
k = exp(-1 * (1/tau) * kt);
k = k / trapz(k);
% deconvolve data
epsilon = 0.1;
PhotometryFields = {'Photometry', 'PhotometryExpFit'};

for pField = PhotometryFields
    deconvSettings = struct('epsilon', [], 'tau', [], 'kernelLength', []);
    for channel = 1:2
        TE.(pField).data(channel).dFdeconv = bpDeconv(TE.(pField).data(channel).dF, k, epsilon, 'none');
        deconvSettings.epsilon(end+1) = epsilon;
        deconvSettings.tau(end+1) = tau;
        deconvSettings.kernelLength(end+1) = kDuration;
    end
    TE.(pField).settings.deconvSettings = deconvSettings;
end