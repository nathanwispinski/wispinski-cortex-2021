%% This script exists to run each of the individual analysis scripts.
% You could, of course, run each of these separately. This is just to
% recreate all of the figures in the paper, as well as the data used for
% statistical analysis in CSV form.

% Reach area
disp("Running reach area analysis...");
analyze_reachArea;

% Accuracy
disp("Running accuracy analysis...");
analyze_accuracy;

% Reaction time
disp("Running reaction time analysis...");
analyze_RT;

disp("Check the folder for the output CSV and PDF files.")