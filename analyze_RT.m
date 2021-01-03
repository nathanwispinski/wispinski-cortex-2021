% analyze_RT.m

% Loads cleanData.mat file generated from cleanData.m
% Analyzes reaction time data time-locked to when the go cue sound was played

% Requires functions:
    % sine5Cost.m
    % sine61Cost.m
    % fminsearchbnd.m
        % https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon
    % Circular Statistics Toolbox for MATLAB
        % https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

% Nathan Wispinski (nathan3@ualberta.ca)
% July 2020

%% Housekeeping
clear all; close all; clc;

% Add fminsearchbnd folder to path
disp("Adding fminsearchbnd folder to path.")
addpath('fminsearchbnd\')

% Add CircStat2012a folder to path
disp("Adding CircStat2012a folder to path.")
addpath('CircStat2012a\')

% Enter location of .mat files saved from cleanData.m script
workingDir = pwd;
cd(workingDir);

% Load cleaned data
load('cleanData.mat');

%% Description of variables from cleanData.mat

% allMT [#of Participants x #of trials]
    % Time from finger lift to screen touch (seconds)
% allRT [#of Participants x #of trials]
    % Time from go cue beep to finger lift (seconds)
% correct [#of Participants x #of trials]
    % Whether participant touched the correct filled-in target
    % (0=incorrect, 1=correct)
% endPos_collapsed [#of Participants x #of trials]
    % Where the cue stopped when participants lifted their finger
    % Converted to 1-61 circle points, where 1=100% , 30=50%, 61=100%
    % Dictated target probability
% trigPos_collapsed [#of Participants x #of trials]
    % Where the cue was when the go beep occured
    % Converted to 1-5 circle trigger points, where 1=100% , 3=50%, 5=100%
% vertOrHorzTargets [#of Participants x #of trials]
    % Whether targets were arranged vertically or horizontally on a trial
    % (0=Vertical, 1=Horizontal)
% reachAreaNorm [#of Participants x #of trials]
    % Single-trial 2D reach area between trajectory and a straight line
    % between start position and average endpoint for the correct target
    % Calculated on 2D plane of interest (vertical or horizontal)
    % Z-scored within participants for each of the potential target
    % directions (i.e., left/right/up/down)
    % Only calculated for correct trials
% rotCond [#of Participants x #of trials]
    % Which rotation condition the cue was on that trial
    % (0=Clockwise, 1=CounterClockwise, 2=Stationary)
% experiment [1 x #of Participants]
    % Which experiment this subject took part in
    % (1=Box, 2=Arrow)
% numSubs_total [double]
    % Total number of subjects used in final analysis
% maxTrials [double]
    % The maximum number of trials a participant could take part in


%% Sine wave model analysis description

% These data are well-described by a sine wave; the cue (box or arrow)
% rotates around a circle, and so the the probability of each target
% fluctuates over time in a sinusoid pattern with fixed frequency 180 deg/s.

% We can fit a sine wave model to single-trial datapoints in each condition for each subject.
% Then the data become 3 fitted sine wave parameters for each condition for each subject.
% We can think of this as a dimensionality reduction technique (e.g., PCA).
% Then we can statistically compare the sine wave parameters that describe these data.

%% Setup sine wave model

% Example initial parameters for sine wave function
% Fixed frequency/period parameter of 180 deg/s
b(1) = 0.75; % sineMean (i.e., vertical shift)
b(2) = 0.2; % sineRange (i.e., amplitude)
b(3) = deg2rad(1); % sineShift (i.e., phase)

% Sine function to fit (5 discrete points from 0pi to 2pi - to line up with data)
sine5Out = @(b) b(1) + (b(2)*sin(b(3) + (linspace(0,2*pi,5)) ));

% Sine function for plotting (120 discrete points from 0pi to 2pi)
sine120Out = @(b) b(1) + (b(2)*sin(b(3) + (linspace(0,2*pi,120)) ));

% Plot example sine wave
figure; plot(1:120,sine120Out(b));


%% Variables for reaction time and corresponding circle position
% Get reaction time for each trial
% Also note corresponding circle position (1:5)
% [2xTrials] matrix for each person for each condition
    % Note: 96 is max # of trials per condition (e.g., vertical and clockwise)

% Initialize matrices for single-trial data
rt_VertStat = nan(numSubs_total,2,96);
rt_VertCW = nan(numSubs_total,2,96);
rt_VertCCW = nan(numSubs_total,2,96);
rt_HorzStat = nan(numSubs_total,2,96);
rt_HorzCW = nan(numSubs_total,2,96);
rt_HorzCCW = nan(numSubs_total,2,96);

for sub = 1:numSubs_total
        % Vertical stationary
        idx = vertOrHorzTargets(sub,:)==0 & rotCond(sub,:)==2 & correct(sub,:)==1;
        rt_VertStat(sub,1,1:sum(idx)) = allRT(sub,idx);
        rt_VertStat(sub,2,1:sum(idx)) = trigPos_collapsed(sub,idx);

        % Vertical clockwise
        idx = vertOrHorzTargets(sub,:)==0 & rotCond(sub,:)==0 & correct(sub,:)==1;
        rt_VertCW(sub,1,1:sum(idx)) = allRT(sub,idx);
        rt_VertCW(sub,2,1:sum(idx)) = trigPos_collapsed(sub,idx);

        % Vertical counterclockwise
        idx = vertOrHorzTargets(sub,:)==0 & rotCond(sub,:)==1 & correct(sub,:)==1;
        rt_VertCCW(sub,1,1:sum(idx)) = allRT(sub,idx);
        rt_VertCCW(sub,2,1:sum(idx)) = trigPos_collapsed(sub,idx);

        % Horizontal stationary
        idx = vertOrHorzTargets(sub,:)==1 & rotCond(sub,:)==2 & correct(sub,:)==1;
        rt_HorzStat(sub,1,1:sum(idx)) = allRT(sub,idx);
        rt_HorzStat(sub,2,1:sum(idx)) = trigPos_collapsed(sub,idx);

        % Horizontal clockwise
        idx = vertOrHorzTargets(sub,:)==1 & rotCond(sub,:)==0 & correct(sub,:)==1;
        rt_HorzCW(sub,1,1:sum(idx)) = allRT(sub,idx);
        rt_HorzCW(sub,2,1:sum(idx)) = trigPos_collapsed(sub,idx);

        % Horizontal counterclockwise
        idx = vertOrHorzTargets(sub,:)==1 & rotCond(sub,:)==1 & correct(sub,:)==1;
        rt_HorzCCW(sub,1,1:sum(idx)) = allRT(sub,idx);
        rt_HorzCCW(sub,2,1:sum(idx)) = trigPos_collapsed(sub,idx);
end


%% Fit sinewave to single-trial datapoints in each condition
% Use fminsearchbnd instead of fminsearch to constrain amplitude parameter to be positive
    % Sometimes instead of shifting 180 degrees, amplitude parameter would switch sign
% NOTE: Parameter fitting should take about 2-5 minutes

% fminsearchbnd options
opts = optimset('fminsearch');
opts.Display = 'off';
LB = [0.100 0.01 -inf]; % Lower bound for each parameter
UB = [0.350 0.200 inf]; % Upper bound for each parameter
nFits = 100; % Number of fits with random parameters to try

% Initialize matrices for fitted parameters
params_VertStat = nan(numSubs_total,3);
params_VertCW = nan(numSubs_total,3);
params_VertCCW = nan(numSubs_total,3);
params_HorzStat = nan(numSubs_total,3);
params_HorzCW = nan(numSubs_total,3);
params_HorzCCW = nan(numSubs_total,3);

for cond = 1:6
    % Display current progress
    disp(['Fitting Condition: ' num2str(cond) '/6']);
for sub = 1:numSubs_total
    
    % Get data for this condition and subject
    if cond==1
        data = squeeze(rt_VertStat(sub,:,:));
    elseif cond==2
        data = squeeze(rt_VertCW(sub,:,:));
    elseif cond==3
        data = squeeze(rt_VertCCW(sub,:,:));
    elseif cond==4
        data = squeeze(rt_HorzStat(sub,:,:));
    elseif cond==5
        data = squeeze(rt_HorzCW(sub,:,:));
    elseif cond==6
        data = squeeze(rt_HorzCCW(sub,:,:));
    end
    
    % Initialize matrix for several fits
    tmp = NaN(nFits,4);
    for iFit = 1:nFits
        % Get random initial parameters within bounds
        randParams = [rand(1) rand(1)*UB(2) rand(1)*deg2rad(360)];
        % Find sine wave parameters by minimizing least-squares to single-trial data
            % Using sine5Cost.m function
        [s, score] = fminsearchbnd(@(x) sine5Cost(x,data), ...
            randParams,LB,UB,opts);
        tmp(iFit,:) = [s score];
    end
    % Get best fitting set of parameters for this subject/condition combo
    bestFit = find(min(tmp(:,4)));
    
    % Save data for this condition and subject
    if cond==1
        params_VertStat(sub,:) = tmp(bestFit,1:3);
    elseif cond==2
        params_VertCW(sub,:) = tmp(bestFit,1:3);
    elseif cond==3
        params_VertCCW(sub,:) = tmp(bestFit,1:3);
    elseif cond==4
        params_HorzStat(sub,:) = tmp(bestFit,1:3);
    elseif cond==5
        params_HorzCW(sub,:) = tmp(bestFit,1:3);
    elseif cond==6
        params_HorzCCW(sub,:) = tmp(bestFit,1:3);
    end
end
end


%% Calculate R2 values for each condition/participant
% Calculate how much variance is accounted for by the sine wave model
% relative to a straight line through the mean
% sine5Cost() outputs the sum of squared deviations between sine model and data

% Initialize r2 matrix
r2_values = nan(numSubs_total,6);

for sub = 1:numSubs_total
    % Horizontal Stationary
    SSModel = sine5Cost(params_HorzStat(sub,:),squeeze(rt_HorzStat(sub,:,:)));
    SSTotal = nansum((squeeze(rt_HorzStat(sub,1,:)) - nanmean(squeeze(rt_HorzStat(sub,1,:)))).^2);
    r2_values(sub,1) = 1-(SSModel/SSTotal);
    
    % Horizontal Clockwise
    SSModel = sine5Cost(params_HorzCW(sub,:),squeeze(rt_HorzCW(sub,:,:)));
    SSTotal = nansum((squeeze(rt_HorzCW(sub,1,:)) - nanmean(squeeze(rt_HorzCW(sub,1,:)))).^2);
    r2_values(sub,2) = 1-(SSModel/SSTotal);
    
    % Horizontal Counterclockwise
    SSModel = sine5Cost(params_HorzCCW(sub,:),squeeze(rt_HorzCCW(sub,:,:)));
    SSTotal = nansum((squeeze(rt_HorzCCW(sub,1,:)) - nanmean(squeeze(rt_HorzCCW(sub,1,:)))).^2);
    r2_values(sub,3) = 1-(SSModel/SSTotal);
    
    % Vertical Stationary
    SSModel = sine5Cost(params_VertStat(sub,:),squeeze(rt_VertStat(sub,:,:)));
    SSTotal = nansum((squeeze(rt_VertStat(sub,1,:)) - nanmean(squeeze(rt_VertStat(sub,1,:)))).^2);
    r2_values(sub,4) = 1-(SSModel/SSTotal);
    
    % Vertical Clockwise
    SSModel = sine5Cost(params_VertCW(sub,:),squeeze(rt_VertCW(sub,:,:)));
    SSTotal = nansum((squeeze(rt_VertCW(sub,1,:)) - nanmean(squeeze(rt_VertCW(sub,1,:)))).^2);
    r2_values(sub,5) = 1-(SSModel/SSTotal);
    
    % Vertical Counterclockwise
    SSModel = sine5Cost(params_VertCCW(sub,:),squeeze(rt_VertCCW(sub,:,:)));
    SSTotal = nansum((squeeze(rt_VertCCW(sub,1,:)) - nanmean(squeeze(rt_VertCCW(sub,1,:)))).^2);
    r2_values(sub,6) = 1-(SSModel/SSTotal);
end


%% Save R2 values to csv (if desired)
csvwrite('R2_RT.csv',r2_values);


%% Convert Phase parameter estimates to [0 360]
% fminsearchbnd Phase bounds are [-inf inf] to prevent estimates converging to bounds
    % We let Phase vary unconstrained, and then convert back to [0 360]
% Find all phase parameters greater than 2pi or less than 0pi
% Wrap around so they fall in [0pi 2pi]

% Horizontal Stationary
idx = params_HorzStat(:,3)>deg2rad(360) | params_HorzStat(:,3)<0;
params_HorzStat(idx,3) = mod(params_HorzStat(idx,3),deg2rad(360));

% Horizontal Clockwise
idx = params_HorzCW(:,3)>deg2rad(360) | params_HorzCW(:,3)<0;
params_HorzCW(idx,3) = mod(params_HorzCW(idx,3),deg2rad(360));

% Horizontal Counterclockwise
idx = params_HorzCCW(:,3)>deg2rad(360) | params_HorzCCW(:,3)<0;
params_HorzCCW(idx,3) = mod(params_HorzCCW(idx,3),deg2rad(360));

% Vertical Stationary
idx = params_VertStat(:,3)>deg2rad(360) | params_VertStat(:,3)<0;
params_VertStat(idx,3) = mod(params_VertStat(idx,3),deg2rad(360));

% Vertical Clockwise
idx = params_VertCW(:,3)>deg2rad(360) | params_VertCW(:,3)<0;
params_VertCW(idx,3) = mod(params_VertCW(idx,3),deg2rad(360));

% Vertical Counterclockwise
idx = params_VertCCW(:,3)>deg2rad(360) | params_VertCCW(:,3)<0;
params_VertCCW(idx,3) = mod(params_VertCCW(idx,3),deg2rad(360));


%% Plot distributions of parameters over 6 conditions
% Add jitter to x coordinate for display purposes
% Stationary = Black
% Clockwise = Blue
% Counterclockwise = Red
% Horizontal = First 3
% Vertical = Last 3

% Make a new figure
figure; hold on;
% Plot Mean parameter estimates
subplot(1,3,1); hold on;
scatter((0.5.*rand(1,numSubs_total))+0.75,params_HorzStat(:,1),'k','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+1.75,params_HorzCW(:,1),'b','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+2.75,params_HorzCCW(:,1),'r','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+4.75,params_VertStat(:,1),'k','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+5.75,params_VertCW(:,1),'b','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+6.75,params_VertCCW(:,1),'r','filled','MarkerEdgeColor','k');
ylim([0 0.350]);
xlim([0.5 7.5]);
title('Mean');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',20);
yticks([0 .150 0.300]); yticklabels({'0' '150' '300'});
xticks([2 6]); xticklabels({'H' 'V'});

% Plot Amplitude parameter estimates
subplot(1,3,2); hold on;
scatter((0.5.*rand(1,numSubs_total))+0.75,params_HorzStat(:,2),'k','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+1.75,params_HorzCW(:,2),'b','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+2.75,params_HorzCCW(:,2),'r','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+4.75,params_VertStat(:,2),'k','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+5.75,params_VertCW(:,2),'b','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+6.75,params_VertCCW(:,2),'r','filled','MarkerEdgeColor','k');
ylim([0 0.100]);
xlim([0.5 7.5]);
title('Amplitude');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',20);
yticks([0 .100 .200]); yticklabels({'0' '100' '200'});
xticks([2 6]); xticklabels({'H' 'V'});

% Plot Phase parameter estimates
subplot(1,3,3); hold on;
scatter((0.5.*rand(1,numSubs_total))+0.75,params_HorzStat(:,3),'k','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+1.75,params_HorzCW(:,3),'b','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+2.75,params_HorzCCW(:,3),'r','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+4.75,params_VertStat(:,3),'k','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+5.75,params_VertCW(:,3),'b','filled','MarkerEdgeColor','k');
scatter((0.5.*rand(1,numSubs_total))+6.75,params_VertCCW(:,3),'r','filled','MarkerEdgeColor','k');
% Plot expected phase
plot([0.5 7.5],[deg2rad(270) deg2rad(270)],'k--','LineWidth',2);
% Plotting options
ylim([deg2rad(0) deg2rad(360)]);
xlim([0.5 7.5]);
title('Phase');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',20);
yticks([0 1*pi 2*pi]); yticklabels({'0' '180' '360'});
xticks([2 6]); xticklabels({'H' 'V'});

% Plotting options
set(gcf,'color','w');


%% Save figure (if desired)
saveas(gcf,'parameters_RT.pdf'); % As .pdf for vector graphics


%% Arrange data for statistics
% For each fitted sine wave parameter seperately
% [numSubs x 6 conditions]

% Mean parameters for both experiments
rt_Mean = ...
[params_HorzStat(:,1) params_HorzCW(:,1) params_HorzCCW(:,1) ...
params_VertStat(:,1) params_VertCW(:,1) params_VertCCW(:,1)];
% Add last column indicating which experiment participant belonged to (Box=1; Arrow=2)
rt_Mean(:,end+1) = experiment;

% Amplitude parameters for both experiments
rt_Amplitude = ...
[params_HorzStat(:,2) params_HorzCW(:,2) params_HorzCCW(:,2) ...
params_VertStat(:,2) params_VertCW(:,2) params_VertCCW(:,2)];
% Add last column indicating which experiment participant belonged to (Box=1; Arrow=2)
rt_Amplitude(:,end+1) = experiment;

% Phase parameters for both experiments
rt_Phase = ...
[params_HorzStat(:,3) params_HorzCW(:,3) params_HorzCCW(:,3) ...
params_VertStat(:,3) params_VertCW(:,3) params_VertCCW(:,3)];
% Add last column indicating which experiment participant belonged to (Box=1; Arrow=2)
rt_Phase(:,end+1) = experiment;


%% Save .csv files of data for statistics elsewhere (if desired)
csvwrite('RT_Mean.csv',rt_Mean);
csvwrite('RT_Amplitude.csv',rt_Amplitude);
csvwrite('RT_Phase.csv',rt_Phase);


%% Circular statistics for comparing phase
% Get grand average phase parameter estimate (collapsing across experiment and arrangement)
% Calculate confidence intervals (w/ Bonferroni correction)
% Test if average phase in each condition is different from expected phase
% NOTE: Requires Circular Statistics Toolbox

% Initialize arrays for test statistics
sigPhase = nan(1,3);
meanPhase = nan(1,3);
upperPhaseCI = nan(1,3);
lowerPhaseCI = nan(1,3);

for cond = 1:3
    % Get phase estimates in this condition
    tmpAngles = circ_mean(rt_Phase(:,[cond cond+3])');
    tmpAngles = wrapTo2Pi(tmpAngles);
    % One-sample test for mean angle against expected angle (w/ confidence intervals)
    [sigPhase(cond), meanPhase(cond), upperPhaseCI(cond), lowerPhaseCI(cond)] = ...
        circ_mtest(tmpAngles,deg2rad(270),0.05/(3*3*2));
end
% Wrap angles to 2pi
meanPhase = wrapTo2Pi(meanPhase);
upperPhaseCI = wrapTo2Pi(upperPhaseCI);
lowerPhaseCI = wrapTo2Pi(lowerPhaseCI);


%% Estimate time lag for each phase

% Estimate mean phase for each rotation collapsing across arrangement
meanPhaseStat = meanPhase(1);
meanPhaseCW = meanPhase(2);
meanPhaseCCW = meanPhase(3);

% Get difference from expected phase
diffPhaseStat = meanPhaseStat - deg2rad(270);
diffPhaseCW = meanPhaseCW - deg2rad(270);
diffPhaseCCW = meanPhaseCCW - deg2rad(270);

% Convert to degrees
diffPhaseStat = rad2deg(diffPhaseStat);
diffPhaseCW = rad2deg(diffPhaseCW);
diffPhaseCCW = rad2deg(diffPhaseCCW);

% Convert degrees to time estimate (stimulus rotates at 180deg/s)
disp(diffPhaseStat/180);
disp(diffPhaseCW/180);
disp(diffPhaseCCW/180);


%% Pairwise comparisons between phase estimates
% Conduct a one-sample circular t-test on the paired differences in each condition

% Initialize arrays for test statistics
sigPhase_paired = nan(3,3);
meanPhase_paired = nan(3,3);
upperPhaseCI_paired = nan(3,3);
lowerPhaseCI_paired = nan(3,3);

for cond = 1:3
    for cond2 = 1:3

    % Get phase estimates in this condition
    tmpAngles1 = wrapTo2Pi(circ_mean([rt_Phase(:,cond) rt_Phase(:,cond+3)]'));
    % Get phase estimates in another condition
    tmpAngles2 = wrapTo2Pi(circ_mean([rt_Phase(:,cond2) rt_Phase(:,cond2+3)]'));
    % Get paired difference of angles
    diffAngles = tmpAngles1 - tmpAngles2;
    % One-sample test for angle difference vs. 0 (w/ confidence intervals)
    [sigPhase_paired(cond,cond2), meanPhase_paired(cond,cond2), ...
        upperPhaseCI_paired(cond,cond2), lowerPhaseCI_paired(cond,cond2)] = ...
        circ_mtest(diffAngles,deg2rad(0),0.05/(3*3*2));
    end
end


%% Make a polar plot
% Average Horizontal and Vertical average phases and upper and lower CI estimates

% Dimensions of the circle to plot
r = 1; % Radius
ang = [0:0.01:2*pi 0]; % Angular points to plot

% Make a new figure
figure; hold on;

% CI patches (corrected for multiple comparisons)
% Stationary
xp = [0 r*cos(linspace(lowerPhaseCI(1)+pi,upperPhaseCI(1)+pi,100)) 0];
yp = [0 r*sin(linspace(lowerPhaseCI(1)+pi,upperPhaseCI(1)+pi,100)) 0];
patch(xp,yp,'k','EdgeColor','None');
% Clockwise
xp = [0 r*cos(linspace(lowerPhaseCI(2)+pi,upperPhaseCI(2)+pi,100)) 0];
yp = [0 r*sin(linspace(lowerPhaseCI(2)+pi,upperPhaseCI(2)+pi,100)) 0];
patch(xp,yp,'b','EdgeColor','None');
% Counterclockwise
xp = [0 r*cos(linspace(lowerPhaseCI(3)+pi,upperPhaseCI(3)+pi,100)) 0];
yp = [0 r*sin(linspace(lowerPhaseCI(3)+pi,upperPhaseCI(3)+pi,100)) 0];
patch(xp,yp,'r','EdgeColor','None');

% Plot mean angles of stationary, CW, and CCW
% Stationary
xp = r*cos(meanPhase(1)+pi);
yp = r*sin(meanPhase(1)+pi);
plot([0 xp],[0 yp],'k','LineWidth',2);
% Clockwise
xp = r*cos(meanPhase(2)+pi);
yp = r*sin(meanPhase(2)+pi);
plot([0 xp],[0 yp],'b','LineWidth',2);
% Counterclockwise
xp = r*cos(meanPhase(3)+pi);
yp = r*sin(meanPhase(3)+pi);
plot([0 xp],[0 yp],'r','LineWidth',2);

% Plot expected angle as dotted line
xp = r*cos(deg2rad(270)+pi);
yp = r*sin(deg2rad(270)+pi);
plot([0 xp],[0 yp+0.2],'k--','LineWidth',2);
% Plot a circle
xp = r*cos(ang);
yp = r*sin(ang);
plot(xp,yp,'k','LineWidth',3);
% Plotting options
alpha(0.3); % Set patch transparency
set(gca,'LineWidth',1.5);
set(gca,'FontSize',20);
axis off; axis equal;
set(gcf,'color','w');


%% Save polar plot of phase parameters
saveas(gcf,'RT_Phase.pdf'); % As .pdf for vector graphics


%% Average parameters to get average sine wave for each Rotation condition

% Collapse parameters across target arrangement conditions (Horizontal/Vertical)
stat_params = [params_HorzStat; params_VertStat];
CW_params = [params_HorzCW; params_VertCW];
CCW_params = [params_HorzCCW; params_VertCCW];

% Average parameters in each rotation condition
avg_statParams = [mean(stat_params(:,1)) mean(stat_params(:,2)) mod(circ_mean(stat_params(:,3)),deg2rad(360))];
avg_CWParams = [mean(CW_params(:,1)) mean(CW_params(:,2)) mod(circ_mean(CW_params(:,3)),deg2rad(360))];
avg_CCWParams = [mean(CCW_params(:,1)) mean(CCW_params(:,2)) mod(circ_mean(CCW_params(:,3)),deg2rad(360))];


%% Plot all individual sine models for each rotation condition
scatterSize = 100;
sineColor = [0.7 0.7 0.7];

% Make a figure
figure; hold on;

% Plot Stationary models
subplot(1,3,1); hold on;
% Plot all individual model fits from stationary conditions
for sub = 1:numSubs_total
    plot(linspace(1,5,120),sine120Out(params_HorzStat(sub,:)),'Color',sineColor);
end
for sub = 1:numSubs_total
    plot(linspace(1,5,120),sine120Out(params_VertStat(sub,:)),'Color',sineColor);
end
% Plot average model
plot(linspace(1,5,120),sine120Out(avg_statParams),'Color','k','LineWidth',2);
% Plotting options
ylim([0.100 0.350]); xlim([0.9 5.1]);
title('Stationary');
xticks([1 3 5]); xticklabels({'100%' '50%' '100%'});
yticks([0.100 0.200 0.300]); yticklabels({'100' '200' '300'});
xlabel('Target Validity');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',20);

% Plot CW data and models
subplot(1,3,2); hold on;
% Plot all individual model fits from stationary conditions
for sub = 1:numSubs_total
    plot(linspace(1,5,120),sine120Out(params_HorzCW(sub,:)),'Color',sineColor);
end
for sub = 1:numSubs_total
    plot(linspace(1,5,120),sine120Out(params_VertCW(sub,:)),'Color',sineColor);
end
% Plot average model
plot(linspace(1,5,120),sine120Out(avg_CWParams),'Color','b','LineWidth',2);
% Plotting options
ylim([0.100 0.350]); xlim([0.9 5.1]);
title('CW');
xticks([1 3 5]); xticklabels({'100%' '50%' '100%'});
yticks([0.100 0.200 0.300]); yticklabels({'100' '200' '300'});
xlabel('Target Validity');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',20);

% Plot CCW data and models
subplot(1,3,3); hold on;
% Plot all individual model fits from stationary conditions
for sub = 1:numSubs_total
    plot(linspace(1,5,120),sine120Out(params_HorzCCW(sub,:)),'Color',sineColor);
end
for sub = 1:numSubs_total
    plot(linspace(1,5,120),sine120Out(params_VertCCW(sub,:)),'Color',sineColor);
end
% Plot average model
plot(linspace(1,5,120),sine120Out(avg_CCWParams),'Color','r','LineWidth',2);
% Plotting options
ylim([0.100 0.350]); xlim([0.9 5.1]);
title('CCW');
xticks([1 3 5]); xticklabels({'100%' '50%' '100%'});
yticks([0.100 0.200 0.300]); yticklabels({'100' '200' '300'});
xlabel('Target Validity');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',20);

% Plotting options
set(gcf,'color','w');


%% Save figure
saveas(gcf,'RT.pdf'); % As .pdf for vector graphics

