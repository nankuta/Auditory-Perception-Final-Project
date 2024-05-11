%% Auditory Perception Final Project Spring 2024
% a dynamical systems model of stream segregation in galloping ABA sequences
% Nick Ankuta
% Prof. David Poeppel

% constructing stimuli from Gutschalk et al. (2005), with semitone differences
% of [2, 4, 6, 8, 10, 12] and interstimulus intervals (ISI's) of 50 vs. 200 ms

% collecting data from 900 (or some number) of subjects

% simulating responses and fitting the supercritical pitchfork bifurcation:
% dy/dx = y(r-y^2)
% where y is the difference between semitones (Δf) and x is the consistency
% with which 2 streams are perceived
% x = 5 reflects 100% consistent perception of 2 streams
% x = 0 reflects 0% consistent perception of 2 streams (or 100% consistent
% perception of 1 stream)

% stable equilibria of the bifurcation are given by y = +/- sqrt(r)
% and an unstable equilibrium exists at y = 0

% goal of the experiment is to recover r

% here the unstable equilibrium corresponds to the Δf point at which
% perceiving 1 vs. 2 streams is equally likely. we can set this point to 0
% with a linear shift. the stable equilibria then correspond to Δf values
% which are +/-sqrt(r) from this point, with -sqrt(r) reflecting the point
% that 1 stream is consistently perceived and +sqrt(r) reflecting the point
% that 2 streams are consistently perceived

% this experiment attempts to formalize categorical perception in the
% galloping illusion and provide a dynamical systems model that predicts
% the consistency of perceiving 1 vs. 2 streams for arbitrary Δf

% this model characterizes the fission / temporal coherence boundary described
% in Van Noorden (1975), and can potentially be applied to future analyses
% of more complex tonal sequences (and answer more elaborate questions
% in categorical perception such as when 1 vs. 2 melodies/voices are perceived
% in a piece of music)

%% first, make the stimuli

% B tone fixed at 1000 Hz
% A tone is 2, 4, 6 ,8, 10, or 12 semitones above the B tone (1122, 1260,
% 1414, 1587, 1782, or 2000 Hz)

% responses are being simulated for this project, but if I were to pursue
% this further, I would construct tonal sequences using the PsychPortAudio
% function of the Psychophysics Toolbox or something similar

%% the experiment

trials_per_condition = 10;
n_subs = 900;

n_semitones = 6;
n_resp_vals = 6;
n_ISI = 2;
n_conditions = n_semitones*n_ISI;

n_trials = trials_per_condition*n_conditions;

semitone_vals = [2, 4, 6, 8, 10, 12];
semitones = zeros(n_semitones, n_ISI, trials_per_condition, n_subs);

for ii = 1:n_semitones
   semitones(ii, :, :, :) = semitone_vals(ii); 
end

% each trial: presentation of a sequence with random semitone difference and
% 50 or 200 ms ISI
% subject asked to rate (0-5) ease of holding 2-stream perception
% 5 reflects 100% consistency, i.e., perceiving only two streams
% 0 reflects 0% consistency, i.e., perceiving only one stream

%% simulate subject responses

% 1) completely random responses
% simulated_resps = randi([0,5], n_semitones, n_ISI, trials_per_condition, n_subs);

% 2) random responses with linear increase in ease of streaming with Δf
simulated_resps = zeros(n_semitones, n_ISI, trials_per_condition, n_subs);
for ii = 0:n_resp_vals-1
   simulated_resps(ii+1, :, :, :) = ii + (1-2*rand(1, n_ISI, trials_per_condition, n_subs));
end

simulated_resps = round(simulated_resps);
simulated_resps = changem(simulated_resps, [0, 5], [-1, 6]);

% My priority for this project was computing r accurately from arbitrary data
% rather than obtaining a good approximation for r itself. I could improve
% my simulated data by using distributions around the experimental values from
% Gutschalk et al. if I continue to pursue this.

%% fit supercritical pitchfork bifurcation

% dy/dx = y(r-y^2);
% bifurcation parameter r characterizes the fission / temporal coherence boundary
% -sqrt(r) reflects the Δf value at which 1 stream is consistently perceived
% and +sqrt(r) reflects the Δf value at which 2 streams are consistently
% perceived
% r can be approximated at a given point using the derivative and the Δf
% value. an average of these approximations is used to estimate the true r value.

% normalize response values
simulated_resps = simulated_resps/5;

% recovered r parameters for each subject for each ISI, initialized to zero
r_params_by_sub = zeros(n_semitones-1, n_ISI, trials_per_condition, n_subs);

% since we are always comparing adjacent points (e.g., Δf = 6 and 8) as we
% compute the derivative, dy is always equal to 2
dy = 2*ones(1, n_ISI, trials_per_condition, n_subs);
y = ones(1, n_ISI, trials_per_condition, n_subs);

% calculate r across trials
for xi = 1:n_semitones-1
    dx = simulated_resps(xi+1, :, :, :) - simulated_resps(xi, :, :, :);
    ydot = dy./dx; % compute derivatives
    y = y+2; % semitone value (midpoint)
    r_params_by_sub(xi, :, :, :) = (ydot + y.^3)./y;
end

% clean data
r_params_by_sub(isinf(r_params_by_sub)) = NaN;

% sort data by ISI
r_params_50 = r_params_by_sub(:, 1, :, :);
r_params_200 = r_params_by_sub(:, 2, :, :);

% recover r
r_param_50_recovered = nanmean(r_params_50, "all");
r_param_200_recovered = nanmean(r_params_200, "all");

% compute standard deviation
r_params_50_std = nanstd(r_params_50, 1, "all");
r_params_200_std = nanstd(r_params_200, 1, "all");

% say that 50-50 point is given by roughly 5 semitones for 50 ms ISI and 7
% semitones for 200 ms ISI
unstable_pt_50 = 5;
unstable_pt_200 = 7;

% then:
crit_1_stream_50 = unstable_pt_50 -sqrt(r_param_50_recovered);
crit_2_stream_50 = unstable_pt_50 + sqrt(r_param_50_recovered);

crit_1_stream_200 = 7 - sqrt(r_param_200_recovered);
crit_2_stream_200 = 7 + sqrt(r_param_200_recovered);
