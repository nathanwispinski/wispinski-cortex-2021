% cleanData.m

% Take .mat files for each participant, clean, and save data for further
% analysis with analyze_*.m scripts.

% Nathan Wispinski (nathan3@ualberta.ca)
% July 2020


%% Housekeeping
clear all; close all; clc;

% Enter location of .mat files with each subject's data
dataDir_Box = 'Data\Box';
dataDir_Arrow = 'Data\Arrow';
% Enter location to save cleaned data
workingDir = pwd;
cd(workingDir);


%% Conditions matrix
% Recreate the conditions matrix used in the experiment script
% This matrix defined each kind of unique trial, and how often it would be
% presented, before being randomized for each participant

% configAndProb = [target1 target2 moveDir moveOnset numTrials]
% For target1 and target2 combinations:
    % 1 and 3 = Horizontal Targets
    % 2 and 4 = Vertical Targets
% Movement Directions:
    % 0 = Clockwise (cw)
    % 1 = CounterClockwise (ccw)
    % 2 = Stationary
% Movement Onset corresponds to each of the 8 trigger locations along the circle
configAndProb = [ ...
1 3 0 1 12; 1 3 1 1 12; 1 3 2 1 12; 2 4 0 1 12; 2 4 1 1 12; 2 4 2 1 12;
1 3 0 2 12; 1 3 1 2 12; 1 3 2 2 12; 2 4 0 2 12; 2 4 1 2 12; 2 4 2 2 12;
1 3 0 3 12; 1 3 1 3 12; 1 3 2 3 12; 2 4 0 3 12; 2 4 1 3 12; 2 4 2 3 12;
1 3 0 4 12; 1 3 1 4 12; 1 3 2 4 12; 2 4 0 4 12; 2 4 1 4 12; 2 4 2 4 12;
1 3 0 5 12; 1 3 1 5 12; 1 3 2 5 12; 2 4 0 5 12; 2 4 1 5 12; 2 4 2 5 12;
1 3 0 6 12; 1 3 1 6 12; 1 3 2 6 12; 2 4 0 6 12; 2 4 1 6 12; 2 4 2 6 12;
1 3 0 7 12; 1 3 1 7 12; 1 3 2 7 12; 2 4 0 7 12; 2 4 1 7 12; 2 4 2 7 12;
1 3 0 8 12; 1 3 1 8 12; 1 3 2 8 12; 2 4 0 8 12; 2 4 1 8 12; 2 4 2 8 12;
];

% Trigger points around circle (for use with 120 xy points around circle)
trigPosFrames = [15 30 45 60 75 90 105 120];


%% Example lines of code from experiment script to help with circle positions
mainRect = [0 0 1920 1080]; % Mock screen dimensions
pxl2HexCent = 243; % 9 cm in pixels (origin to target center)
pxlTargetRadius = 27; % 1 cm in pixels (radius of target)
a = mainRect(3)/2; % Center of screen in x pixels
b = mainRect(4)/2; % Center of screen in y pixels
r = pxl2HexCent-2*pxlTargetRadius; % Draw box so far edge just touches target (7 cm from origin)
t = linspace(0,2*pi,120); % 120 frames at 60 Hz means 500ms per quadrant
% Mock x and y coordinates
x = a + r*cos(t);
y = b + r*sin(t);

% OPTIONAL: Plot all 120 unique circle positions for this mock setup
figure; scatter(x,y); axis equal;


%% Initialize variables to put participant data into

% Make sure we're in working directory
cd(workingDir);

% Get names of all files in Box experiment
cd(dataDir_Box);
subFiles_Box = dir('*.mat'); % Find all .mat files in current folder
numSubs_Box = length(subFiles_Box); % How many subjects are there
cd(workingDir); % Back to working directory

% Get names of all files in Arrow experiment
cd(dataDir_Arrow);
subFiles_Arrow = dir('*.mat'); % Find all .mat files in current folder
numSubs_Arrow = length(subFiles_Arrow); % How many subjects are there
cd(workingDir); % Back to working directory

% Get total number of usable subjects in both experiments
numSubs_total = numSubs_Box + numSubs_Arrow;

% Maximum number of possible trials for each participant
maxTrials = 576;

% Initialize some variables
% Get condition #s for each trial
trialConditions = nan(numSubs_total,maxTrials);
% Find whether targets on each trial were vertical or horizontal
vertOrHorzTargets = nan(numSubs_total,maxTrials);
% Get trigger position of each trial
trigPos = nan(numSubs_total,maxTrials);
% Find if/which way box was moving
rotCond = nan(numSubs_total,maxTrials);
% Get accuracy for each trial
correct = nan(numSubs_total, maxTrials);
% Get all RTs for each trial
allRT = nan(numSubs_total, maxTrials);
% Get all MTs for each trial
allMT = nan(numSubs_total, maxTrials);
% Get final box position for each trial (numTrials x,y);
finalBoxPos = nan(numSubs_total,maxTrials,2);
% Get 3D reach trajectories
xReaches = nan(numSubs_total, maxTrials, 200);
yReaches = nan(numSubs_total, maxTrials, 200);
zReaches = nan(numSubs_total, maxTrials, 200);
% Get trials where the only error was tooSlow
idxTooSlow = nan(numSubs_total, maxTrials);
% Get counts for errors
count_numTrials = nan(1,numSubs_total);
count_tooEarly = nan(1,numSubs_total);
count_miss = nan(1,numSubs_total);
count_timeOut = nan(1,numSubs_total);


%% Extract data from all datasets for this experiment

% Loop through all subject files
for sub = 1:numSubs_total
    
    % Load Box participants first, then Arrow participants
    if sub <= numSubs_Box
        % Load subject's .mat data file
        cd(dataDir_Box);
        disp(['Box Experiment. Loading ' subFiles_Box(sub).name]);
        load(subFiles_Box(sub).name);
    elseif sub > numSubs_Box
        % Load subject's .mat data file
        cd(dataDir_Arrow);
        disp(['Arrow Experiment. Loading ' subFiles_Arrow(sub-numSubs_Box).name]);
        load(subFiles_Arrow(sub-numSubs_Box).name);
    end
    
    % Find how many trials for this participant (max 576)
    numTrials = length(data.rxnTime);
    
    % Find the condition on each trial
    for i = 1:numTrials
        trialConditions(sub,i) = data.preConditions{i};
    end
    
    % Target arrangement (0=vertical, 1=horizontal)
    vertOrHorzTargets(sub,1:numTrials) = configAndProb(trialConditions(sub,1:numTrials),1) == 1;
    
    % Trigger position
    trigPos(sub,1:numTrials) = configAndProb(trialConditions(sub,1:numTrials),4);
    
    % Rotation condition (0=cw, 1=ccw, 2=stationary)
    rotCond(sub,1:numTrials) = configAndProb(trialConditions(sub,1:numTrials),3);

    % Index correct trials
    correct(sub,1:numTrials) = all(data.error'==0);
    
    % Get all RTs, MTs for each participant
    allRT(sub,1:numTrials) = data.rxnTime;
    allMT(sub,1:numTrials) = data.mvmtTime;
    
    % Get final box position for every trial
    for i = 1:numTrials
        finalBoxPos(sub,i,:) = [data.finalBoxPos{i}(1) data.finalBoxPos{i}(2)];
    end
    
    % Get reach trajectories
    xReaches(sub,1:size(data.fdaMat.x,1),:) = data.fdaMat.x;
    yReaches(sub,1:size(data.fdaMat.x,1),:) = data.fdaMat.y;
    zReaches(sub,1:size(data.fdaMat.x,1),:) = data.fdaMat.z;
    
    % Errors [1x4]
    % [tooEarly tooSlow Miss timeOut]
    % For accuracy analysis, remove trials where tooEarly and timeOut (i.e., RT too short or too long)
    % Keeping tooEarly and timeOut trials artificially makes accuracy in some
    % probability bins 0, because probability is dictated by RT
    % A side effect of this is that accuracy now only reflects if participants 
    % sucessfully touched the correct target before the movement time limit
    % on trials when participants initiated their reach within RT time limits
    correct(sub,any(data.error(:,[1 4]),2)) = NaN;
    
    % Index tooSlow trials (we will reject or keep them later)
    idxTooSlow(sub,1:size(data.error,1)) = ismember(data.error,[0 1 0 0],'rows');
    
    % Count number of trials, and count each error
    count_numTrials(sub) = numTrials;
    count_tooEarly(sub) = nansum(data.error(:,1));
    count_miss(sub) = nansum(data.error(:,3));
    count_timeOut(sub) = nansum(data.error(:,4));

    % Clear this subject's data and move on to next subject
    clear data
    close all;
    cd(workingDir); % Back to working directory
end
disp('Done.');
cd(workingDir); % Back to working directory


%% Remove trials with 2SD above MT mean

count_badMT = zeros(1,numSubs_total);
for sub = 1:numSubs_total
    
    % Make all tooSlow trials correct
    correct(sub,idxTooSlow(sub,:)==1) = 1;
    
    % Get rid of any trials where MT > 850 ms
    idxBadMT = allMT(sub,:) > 0.850;
    correct(sub,idxBadMT) = 0;
    
    % Then calculate 2SD threshold for MT
    meanMT = nanmean(allMT(sub,correct(sub,:)==1));
    stdMT = nanstd(allMT(sub,correct(sub,:)==1));
    threshMT = meanMT + (2*stdMT);
    
    % Mark all trials as above new MT threshold as incorrect
    % AKA these are the new TooSlow trials
    idxBadMT = allMT(sub,:) > threshMT;
    correct(sub,idxBadMT) = 0;
    
    % Count number of trials that had bad MTs
    count_badMT(sub) = sum(idxBadMT);
end


%% Remove all trials with NaNs in reaches
% NaNs indicate unexpected results from reach cleaning and normalization

count_nanReaches = zeros(1,numSubs_total);
for sub = 1:numSubs_total
    for trial = 1:maxTrials
        % Find if there are any NaNs in x, y, or z for this trial
        xNan = any(isnan(squeeze(xReaches(sub,trial,:))));
        yNan = any(isnan(squeeze(yReaches(sub,trial,:))));
        zNan = any(isnan(squeeze(zReaches(sub,trial,:))));
        
        % If any NaNs, remove trial from analysis
        % Trials that are neither incorrect (0) or correct(0) are not
        % included in analyze_*.m script analysis
        if xNan || yNan || zNan
            correct(sub,trial) = NaN;
            
            % Count number of trials that had bad reaches
            count_nanReaches(sub) = count_nanReaches(sub) + 1;
        end
    end
end


%% Remove trials with bad start or end positions
% NOTE: Reaches are in mm

count_badEndPos = zeros(1,numSubs_total);
for sub = 1:numSubs_total
    
    % Get 3D distance from first reach point to origin [0,0,0]
    startDist = ...
        sqrt(squeeze(xReaches(sub,:,1)).^2 + squeeze(yReaches(sub,:,1)).^2 ...
        + squeeze(zReaches(sub,:,1)).^2);
    % Remove any trials where first reach point is >2cm from origin
    idxBadStartPos = startDist > 20;
    correct(sub,idxBadStartPos) = NaN;
    
    % Get 1D distance at end of reach on towards/away from screen dimension
    endDist = squeeze(yReaches(sub,:,end));
    % Remove any trials where reach endpoint is <37 cm or >43 cm in y dimension
    idxBadEndPos = endDist < 370 | endDist > 430;
    correct(sub,idxBadEndPos) = NaN;
    
    % Count number of trials that had bad start/end positions
    count_badEndPos(sub) = sum(idxBadEndPos);
end


%% Calculate left/right up/down choices based on reach trajectory endpoints

% Initialize matrix for side chosen data
    % On horizontal trials, 1=Left, 2=Right
    % On vertical trials, 1=Down, 2=Up
sideChosen = nan(numSubs_total,maxTrials);
for sub = 1:numSubs_total
    idx = vertOrHorzTargets(sub,:)==1; % Horizontal trials
    sideChosen(sub,idx) = (xReaches(sub,idx,end)>0)+1; % left=1, right=2
    idx = vertOrHorzTargets(sub,:)==0; % Vertical trials
    sideChosen(sub,idx) = (zReaches(sub,idx,end)>250)+1; % down=1, up=2
end


%% Calculate area under curve for each correct reach
% For each subject, calculate average correct endpoint for each target (up/down/left/right)
% For each correct reach, calculate 2D area between reach and "optimal" reach
% where "optimal" reach is a straight line between the start point [0,0,0]
% and the mean correct endpoint for that trial
% Calculate 2D area on plane of interest (e.g., on Horizontal trials,
% calculate area along left/right plane)

% Initialize matrices
reachArea = nan(numSubs_total,numTrials);

for sub = 1:numSubs_total
    
    % Calculate this subject's mean end point for up/down/left/right
    endPos = nan(4,3);
    % Down
    idx = correct(sub,:)==1 & sideChosen(sub,:)==1 & vertOrHorzTargets(sub,:)==0;
    endPos(1,:) = [nanmean(xReaches(sub,idx,end))  nanmean(yReaches(sub,idx,end))  nanmean(zReaches(sub,idx,end))];
    % Up
    idx = correct(sub,:)==1 & sideChosen(sub,:)==2 & vertOrHorzTargets(sub,:)==0;
    endPos(2,:) = [nanmean(xReaches(sub,idx,end))  nanmean(yReaches(sub,idx,end))  nanmean(zReaches(sub,idx,end))];
    % Left
    idx = correct(sub,:)==1 & sideChosen(sub,:)==1 & vertOrHorzTargets(sub,:)==1;
    endPos(3,:) = [nanmean(xReaches(sub,idx,end))  nanmean(yReaches(sub,idx,end))  nanmean(zReaches(sub,idx,end))];
    % Right
    idx = correct(sub,:)==1 & sideChosen(sub,:)==2 & vertOrHorzTargets(sub,:)==1;
    endPos(4,:) = [nanmean(xReaches(sub,idx,end))  nanmean(yReaches(sub,idx,end))  nanmean(zReaches(sub,idx,end))];

    % Get reach area for each trial
    for i = 1:numTrials
        % Get end position for this trial
        if sideChosen(sub,i)==1 && vertOrHorzTargets(sub,i)==0
            curEndPos = endPos(1,:);
        elseif sideChosen(sub,i)==2 && vertOrHorzTargets(sub,i)==0
            curEndPos = endPos(2,:);
        elseif sideChosen(sub,i)==1 && vertOrHorzTargets(sub,i)==1
            curEndPos = endPos(3,:);
        elseif sideChosen(sub,i)==2 && vertOrHorzTargets(sub,i)==1
            curEndPos = endPos(4,:);
        end
        
        % Only for correct trials
        if correct(sub,i)==1
            
            % Get each of the y points for this reach
            yPoints = squeeze(yReaches(sub,i,:));
            
            % If vertical calculate vertical area
            if vertOrHorzTargets(sub,i)==0
                % Make corresponding y points for the "optimal" reach
                slope = curEndPos(3)/curEndPos(2);
                optimalYPoints = slope .* yPoints;
                % Subtract "optimal" reach from actual reach
                flattenedCurve = squeeze(zReaches(sub,i,:)) - optimalYPoints;
            % If horizontal, calculate lateral area
            elseif vertOrHorzTargets(sub,i)==1
                % Make corresponding y points for the "optimal" reach
                slope = curEndPos(1)/curEndPos(2);
                optimalYPoints = slope .* yPoints;
                % Subtract "optimal" reach from actual reach
                flattenedCurve = squeeze(xReaches(sub,i,:)) - optimalYPoints;
            end
            
            % Reverse sign of area on opposite side of space so that
            % positive area always means reaches deviated between the two targets
            if sideChosen(sub,i)==2
                flattenedCurve = flattenedCurve .* -1;
            end
            % Approximate area under curve
            reachArea(sub,i) = trapz(flattenedCurve);
        end
    end
end


%% Normalize reach area for each of the 4 end positions
% Because left/right/up/down reach areas all have different baselines due
% to biomechanics (i.e., right-handed participants reach for left and right
% targets differently)
% Z-score reach areas within each end positions (left/right/up/down) within
% each participant

% Initialize space for normalized reach areas
reachAreaNorm = NaN(size(reachArea));

for sub = 1:numSubs_total
    % Get trial numbers for each reach direction
    idxDown = sideChosen(sub,:)==1 & vertOrHorzTargets(sub,:)==0;
    idxUp = sideChosen(sub,:)==2 & vertOrHorzTargets(sub,:)==0;
    idxLeft = sideChosen(sub,:)==1 & vertOrHorzTargets(sub,:)==1;
    idxRight = sideChosen(sub,:)==2 & vertOrHorzTargets(sub,:)==1;
    
    % Z-score reach areas within each reach direction for this participant
    reachAreaNorm(sub,idxDown) = ...
        (reachArea(sub,idxDown) - nanmean(reachArea(sub,idxDown))) ./ nanstd(reachArea(sub,idxDown));
    reachAreaNorm(sub,idxUp) = ...
        (reachArea(sub,idxUp) - nanmean(reachArea(sub,idxUp))) ./ nanstd(reachArea(sub,idxUp));
    reachAreaNorm(sub,idxLeft) = ...
        (reachArea(sub,idxLeft) - nanmean(reachArea(sub,idxLeft))) ./ nanstd(reachArea(sub,idxLeft));
    reachAreaNorm(sub,idxRight) = ...
        (reachArea(sub,idxRight) - nanmean(reachArea(sub,idxRight))) ./ nanstd(reachArea(sub,idxRight));
end


%% Get circle position when cue ended (of 120 unique positions)

% Initialze space for end position of cue
endPos = nan(numSubs_total,maxTrials);

for sub = 1:numSubs_total
    % 120 unique circle positions for this subject [x,y]
    tmp = squeeze(finalBoxPos(sub,:,:));
    % Estimate 120 unique circle positions for this subject using their circle positions
    radius = mean([max(tmp(:,1))-min(tmp(:,1)); max(tmp(:,2))-min(tmp(:,2))]) /2;
    xCent = (max(tmp(:,1))-min(tmp(:,1)))/2;
    xCent = xCent + min(tmp(:,1));
    yCent = (max(tmp(:,2))-min(tmp(:,2)))/2;
    yCent = yCent + min(tmp(:,2));
    t = linspace(0,2*pi,120);
    x = xCent + radius*cos(t);
    y = yCent + radius*sin(t);
    circlePos = [x;y]';

    % Assign circle positions (1-120) to individual trials
        % Note: circle positions don't overlap perfectly
    for i = 1:size(circlePos,1)
        % Get x and y distance to this circle
        posDev = abs(tmp - circlePos(i,:));
        % Find all trials that are close enough to this circle point
        idx = all(posDev<0.1,2);
        % Assign circle numbers to each trial
        endPos(sub,idx) = i;
    end
end


%% Convert circle positions to probabilies

% Initialize space for circle probabilities
circleProbs = nan(size(circlePos));

% Calculate target probability for each circle point on horizontal trials
circleProbs(:,1) = (circlePos(:,1)-min(circlePos(:,1))) / (max(circlePos(:,1))-min(circlePos(:,1)));
% Calculate target probability for each circle point on vertical trials
circleProbs(:,2) = (circlePos(:,2)-min(circlePos(:,2))) / (max(circlePos(:,2))-min(circlePos(:,2)));


%% Rotate vertical data to align with horizontal data

% For circle positions (1-120)
for sub = 1:numSubs_total
    % Get all vertical trials
    idx = vertOrHorzTargets(sub,:)==0;
    % Subtract 90 degrees (i.e., 30 circle points)
    endPos(sub,idx) = endPos(sub,idx)-30;
    % Get all non-positive end position points
    idx = endPos(sub,:)<1;
    % Rotate so that non-positive circle points wrap back around (i.e., -29 becomes 91)
    endPos(sub,idx) = endPos(sub,idx) + 120;
end

% For trigger positions (1-8)
for sub = 1:numSubs_total
    % Get all vertical trials
    idx = vertOrHorzTargets(sub,:)==0;
    % Subtract 90 degrees (i.e., 2 trigger position points)
    trigPos(sub,idx) = trigPos(sub,idx)-2;
    % Get all non-positive trigger points
    idx = trigPos(sub,:)<1;
    % Rotate so that non-positive circle points wrap back around (i.e., -1 becomes 8)
    trigPos(sub,idx) = trigPos(sub,idx) + 8;
end


%% Collapse endPos into Probability (100% - 50% - 100%)

% For circle positions (1-120)
% Rotate circle points back 1 step, so targets are at 1 and 61
for sub = 1:numSubs_total
    % Add 1 circle position
    endPos(sub,:) = endPos(sub,:)+1;
    % Get all circle points greater than 120 (i.e., max)
    idx = endPos(sub,:)>120;
    % Rotate so that last circle point becomes first circle point
    endPos(sub,idx) = 1;
end
% Take circle points 62-120 and wrap back around to 2-60
endPos_collapsed = NaN(size(endPos));
for sub = 1:numSubs_total
    idx = ismember(endPos(sub,:),62:120);
    endPos_collapsed(sub,:) = endPos(sub,:);
    endPos_collapsed(sub,idx) = endPos_collapsed(sub,idx)-60;
end


% For trigger positions (1-8)
% Rotate circle points back 1 step, so targets are at 1 and 5
for sub = 1:numSubs_total
    % Add 1 trigger position
    trigPos(sub,:) = trigPos(sub,:)+1;
    % Get all trigger points greater than 8 (i.e., max)
    idx = trigPos(sub,:)>8;
    % Rotate so that last trigger point becomes first trigger point
    trigPos(sub,idx) = 1;
end
% Collapse circle points into 1-5 (i.e., target to target)
    % Make trigger points 6 7 8 into 2 3 4
trigPos_collapsed = NaN(size(trigPos));
for sub = 1:numSubs_total
    idx = ismember(trigPos(sub,:),6:8);
    trigPos_collapsed(sub,:) = trigPos(sub,:);
    trigPos_collapsed(sub,idx) = trigPos_collapsed(sub,idx)-4;
end


%% Make variable identifying which participants belong to which experiment
% 1=Box, 2=Arrow
experiment = [ones(numSubs_Box,1); 1+ones(numSubs_Arrow,1)];


%% Save all data for use in analyze_*.m scripts

% Into directory to save .mat file with cleaned data
cd(workingDir);

% Save all relevant variables as cleanData.mat
save('cleanData.mat',...
    'numSubs_total',...
    'maxTrials',...
    'experiment',...
    'vertOrHorzTargets',...
    'rotCond',...
    'correct',...
    'reachAreaNorm',...
    'endPos_collapsed',...
    'trigPos_collapsed',...
    'allMT',...
    'allRT');
disp('Cleaned data saved.');

