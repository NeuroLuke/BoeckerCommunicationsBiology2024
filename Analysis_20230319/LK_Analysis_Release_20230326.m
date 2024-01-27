%==========================================================================
% This script analyzes the behavioral data.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; clc; close all;

% paths
paths       = [];
paths.beh   = 'G:\My Drive\Manuscripts\2023_Boecker_SportsAndSpatialMemory\Preprocessing_20230319\';
paths.save  = 'G:\My Drive\Manuscripts\2023_Boecker_SportsAndSpatialMemory\Analysis_20230319\';

% settings
param                   = [];
param.maxNumTrials      = 160;
param.numTrialsPerChunk = 20;

% logfiles
subjects    = dir(strcat(paths.beh, '*-*'));
fprintf('Number of subjects/sessions: %d.\n', size(subjects, 1));

%% random arena locations for converting the drop errors into memory scores

% possible random locations within arena
dt          = [];
dt.maxR     = 5000; % maximum radius
dt.minR     = 0; % minimum radius
dt.N        = 1000001; % number of locations to create
dt.centerX  = 0; % arena center x
dt.centerY  = 0; % arena center y
locsInArena = LK_RandomPointsInCircle(dt);

% preallocate main output
allRes      = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)

    %% data

    % report
    fprintf('\nAnalysis of session: %s.\n', subjects(iSub).name);

    % load trial data
    trials  = load(fullfile(subjects(iSub).folder, subjects(iSub).name, 'trialInfo.mat'));
    trials  = trials.trialInfo;

    % restrict to maximum number of trials
    fprintf('Number of trials: %d.\n', size(trials, 1));
    if size(trials, 1) > param.maxNumTrials
        trials  = trials(1:param.maxNumTrials, :);
    end
    fprintf('Number of trials after restriction: %d.\n', size(trials, 1));
    numTrials   = size(trials, 1);

    % load behavioral data
    beh = load(fullfile(subjects(iSub).folder, subjects(iSub).name, 'behInfo.mat'));
    beh = beh.behInfo;

    %% task duration

    % duration of the task
    taskStart       = min(trials{1, {'ITI', 'Cue', 'Retrieval', 'Feedback', 'Reencoding', 'Grab'}});
    taskEnd         = max(trials{end, {'ITI', 'Cue', 'Retrieval', 'Feedback', 'Reencoding', 'Grab'}});
    taskDuration    = taskEnd - taskStart;
    fprintf('Task duration: %.3f s.\n', taskDuration);

    %% drop error and memory performance per trial

    % mean drop error, computed by Unreal
    meanDropErrorUnreal     = mean(trials.DropError, 'omitnan');

    % drop error and memory score per trial
    dropError           = sqrt((trials.xCorrect - trials.xResponse) .^ 2 + (trials.yCorrect - trials.yResponse) .^ 2);
    dropErrorSurro      = pdist2([trials.xCorrect, trials.yCorrect], locsInArena);
    MS                  = sum(dropError < dropErrorSurro, 2) / size(dropErrorSurro, 2);

    %% mean drop error and mean memory performance

    % average across all trials
    meanDropError   = mean(dropError, 'omitnan');
    meanMS          = mean(MS, 'omitnan');
    fprintf('Mean drop error: %.3f, mean memory score: %.3f.\n', meanDropError, meanMS);

    %% drop error and memory performance in chunks of trials

    % mean drop error and mean memory score in chunks of data
    chunkedDropError    = nan(ceil(size(dropError, 1) / param.numTrialsPerChunk), 1);
    chunkedMS           = nan(size(chunkedDropError));
    for iChunk = 1:size(chunkedDropError, 1)
        idx                         = [((iChunk - 1) * param.numTrialsPerChunk + 1), (iChunk * param.numTrialsPerChunk)];
        bThisChunk                  = trials.TrialIdx >= min(idx) & trials.TrialIdx <= max(idx);
        chunkedDropError(iChunk, 1) = mean(dropError(bThisChunk), 'omitnan');
        chunkedMS(iChunk, 1)        = mean(MS(bThisChunk), 'omitnan');
    end

    %% retrieval distance and retrieval duration

    % retrieval duration for each trial
    retDur      = trials.Feedback - trials.Retrieval;
    meanRetDur  = mean(retDur, 'omitnan');

    % retrieval path distance for each trial
    retPathDist = nan(size(trials, 1), 1);
    for iTrial = 1:size(trials, 1)

        % behavior during the retrieval period of this trial
        bThisTrial      = beh.trialIdx == iTrial & beh.trialPhase == 3;
        behThisTrial    = beh(bThisTrial, :);

        % retrieval path distance
        retPathDist(iTrial, 1)  = sum(behThisTrial.distances);
    end
    meanRetPathDist = mean(retPathDist, 'omitnan');

    %% correlation between object distance to center and drop error or memory score

    % object distance to the arena center
    objD2Ctr                = sqrt(trials.xCorrect .^ 2 + trials.yCorrect .^ 2);

    % correlation between drop error or memory score and object distance to
    % center
    rhoDropErrorVSObjD2Ctr  = corr(objD2Ctr, dropError, 'type', 'spearman');
    rhoMSVSObjD2Ctr         = corr(objD2Ctr, MS, 'type', 'spearman');

    %% learning rate according to power series

    % fit power series to data
    [fitobject, gof]    = fit(trials.Feedback, MS, 'power1'); % fit uses this function for fitting: y = a * x ^ b

    % learning rate
    learnRate           = fitobject.b;
    learnRateR2         = gof.adjrsquare;

    %% subject distance to the environment center

    % for each time point, estimate the subject's distance to the
    % environment center
    playerD2Ctr     = sqrt(beh.x .^ 2 + beh.y .^ 2);
    
    % mean player distance to center during all cue, retrieval, feedback,
    % and reencoding periods
    bOI                     = beh.trialPhase > 1 & beh.trialIdx <= param.maxNumTrials;
    meanPlayerD2Ctr         = mean(playerD2Ctr(bOI));

    % mean player distance to center during all retrieval periods
    bOI                     = beh.trialPhase == 3 & beh.trialIdx <= param.maxNumTrials;
    meanPlayerD2CtrRet      = mean(playerD2Ctr(bOI));

    % mean player distance to center during all retrieval periods while
    % moving
    bOI                     = beh.trialPhase == 3 & beh.trialIdx < param.maxNumTrials & beh.distances > 0;
    meanPlayerD2CtrRetMove  = mean(playerD2Ctr(bOI));

    %% collect results

    % results from this sessions
    thisRes                         = [];
    thisRes.subject                 = subjects(iSub).name;
    thisRes.numTrials               = numTrials;
    thisRes.meanDropErrorUnreal     = meanDropErrorUnreal;
    thisRes.meanDropError           = meanDropError;
    thisRes.meanMS                  = meanMS;
    thisRes.chunkedDropError        = chunkedDropError;
    thisRes.chunkedMS               = chunkedMS;
    thisRes.taskDuration            = taskDuration;
    thisRes.meanRetDur              = meanRetDur;
    thisRes.meanRetPathDist         = meanRetPathDist;
    thisRes.rhoDropErrorVSObjD2Ctr  = rhoDropErrorVSObjD2Ctr;
    thisRes.rhoMSVSObjD2Ctr         = rhoMSVSObjD2Ctr;
    thisRes.learnRate               = learnRate;
    thisRes.learnRateR2             = learnRateR2;
    thisRes.meanPlayerD2Ctr         = meanPlayerD2Ctr;
    thisRes.meanPlayerD2CtrRet      = meanPlayerD2CtrRet;
    thisRes.meanPlayerD2CtrRetMove  = meanPlayerD2CtrRetMove;

    % collect results across sessions
    allRes                          = cat(1, allRes, thisRes);
end

%% save results

% as matlab file
save(strcat(paths.save, 'allRes'), 'allRes');

% as excel file
data4Excel  = struct2table(allRes);
writetable(data4Excel, strcat(paths.save, 'LK_SportsAndSpatialMemory_20230326.xlsx'));

%% number of trials

% number of trials per session
allNumTrials    = [allRes.numTrials]';
fprintf('\nNumber of trials across sessions: min = %d, max = %d.\n', min(allNumTrials), max(allNumTrials));

%% relationship between mean drop error and mean memory score

% mean drop error and mean memory score for all subjects * sessions
allMeanDropError    = [allRes.meanDropError]';
allMeanMS           = [allRes.meanMS]';

% create figure
f = figure('units', 'centimeters', 'Position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'Position', [1.7, 1.5, 4, 4]);
plot(allMeanDropError, allMeanMS, 'x');
xl = xlabel('Mean drop error (v.u.)');
yl = ylabel('Mean memory score');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
print(f, strcat(paths.save, 'relationMeanDropError2MeanMS'), '-dpng', '-r300');

%% test whether the memory scores are above chance

% test whether the memory scores are above 0.5
[~, p, ~, stats] = ttest(allMeanMS, 0.5);
fprintf('\nComparison of the memory scores against 0.5 chance level: t(%d) = %.3f, p = %.3f.\n', stats.df, stats.tstat, p);

%% chunked drop error and chunked memory score

% all chunked drop errors and memory scores
allChunkedDropError     = transpose(cell2mat({allRes.chunkedDropError}));
allChunkedMS            = transpose(cell2mat({allRes.chunkedMS}));

% plot the average drop error and memory score per chunk
groups  = {'allChunkedDropError', 'Mean drop error'; ...
    'allChunkedMS', 'Mean memory score'};
for iG = 1:size(groups, 1)

    % select data
    if strcmp(groups{iG, 1}, 'allChunkedDropError')
        data    = allChunkedDropError;
    elseif strcmp(groups{iG, 1}, 'allChunkedMS')
        data    = allChunkedMS;
    end

    % create figure
    f = figure('units', 'centimeters', 'Position', [2, 2, 6, 6]);
    axes('units', 'centimeters', 'Position', [1.7, 1.5, 4, 4]);
    bar(mean(data, 1, 'omitnan'), 'FaceColor', [1, 1, 1]);
    xl = xlabel('Chunk');
    yl = ylabel(groups{iG, 2});
    set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
    set(gca, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
    print(f, strcat(paths.save, groups{iG, 1}), '-dpng', '-r300');
end

%% task duration, retrieval duration, and retrieval distance

% task duration
allTaskDur          = [allRes.taskDuration]';
fprintf('\nTask duration: min = %.3f, max = %.3f s.\n', min(allTaskDur), max(allTaskDur));

% retrieval duration
allMeanRetDur       = [allRes.meanRetDur]';
fprintf('Mean retrieval duration: min = %.3f, max = %.3f s.\n', min(allMeanRetDur), max(allMeanRetDur));

% retrieval path distance
allMeanRetPathDist  = [allRes.meanRetPathDist]';
fprintf('Mean retrieval path distance: min = %.3f, max = %.3f s.\n', min(allMeanRetPathDist), max(allMeanRetPathDist));

%% relationship between drop error or memory score and object distance to arena center

% all relationships between drop error or memory score and object distance
% to arena center
allRhoDropErrorVSObjD2Ctr   = [allRes.rhoDropErrorVSObjD2Ctr]';
allRhoMSVSObjD2Ctr          = [allRes.rhoMSVSObjD2Ctr]';

% create figure
f = figure('units', 'centimeters', 'Position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'Position', [1.3, 1.5, 4, 4]);
histogram(allRhoMSVSObjD2Ctr, 'facecolor', [1, 1, 1]);
xl = xlabel('MS vs. object-center distance');
yl = ylabel('Count');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
print(f, strcat(paths.save, 'allRhoMSVSObjD2Ctr'), '-dpng', '-r300');

%% learning rate

% learning rates and their goodness of fit (adjusted r squared)
allLearnRate    = [allRes.learnRate]';
allLearnRateR2  = [allRes.learnRateR2]';

% histogram: all learning rates
f = figure('units', 'centimeters', 'Position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'Position', [1.5, 1.5, 4, 4]);
histogram(allLearnRate, 'facecolor', [1, 1, 1]);
xl = xlabel('Learning rate');
yl = ylabel('Count');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
print(f, strcat(paths.save, 'allLearnRate'), '-dpng', '-r300');

%% player distance to center

% player distances to the center for different conditions
allMeanPlayerD2Ctr          = [allRes.meanPlayerD2Ctr]';
allMeanPlayerD2CtrRet       = [allRes.meanPlayerD2CtrRet]';
allMeanPlayerD2CtrRetMove   = [allRes.meanPlayerD2CtrRetMove]';

% histogram across sessions
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'Position', [1.5, 1.5, 4, 4]);
histogram(allMeanPlayerD2Ctr, 1400:25:3200, 'FaceColor', 'none');
xl = xlabel('Mean distance to center');
yl = ylabel('Count');
set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.4);
set(gca, 'tickdir', 'out', 'ticklength', [0.02, 0.02], 'box', 'off');
print(f, strcat(paths.save, 'allMeanPlayerD2Ctr'), '-dpng', '-r300');

%% reliability

% unique subjects and visits
allSubjects         = {allRes.subject}';
allUniqueSubjects   = unique(cellfun(@(x) x(1), cellfun(@(x) split(x, '-'), allSubjects, 'UniformOutput', false)));
allUniqueVisits     = unique(cellfun(@(x) x(2), cellfun(@(x) split(x, '-'), allSubjects, 'UniformOutput', false)));

% for each unique subject and visit, identify the mean drop error
allMeanDropErrorPerVisit    = nan(size(allUniqueSubjects, 1), size(allUniqueVisits, 1));
for iUS = 1:size(allUniqueSubjects, 1)
    for iUV = 1:size(allUniqueVisits, 1)
        idx = contains(allSubjects, allUniqueSubjects{iUS}) & contains(allSubjects, ['-', allUniqueVisits{iUV}, '-']);
        if sum(idx) > 0
            allMeanDropErrorPerVisit(iUS, iUV) = allMeanDropError(idx);
            % report
            fprintf('Extracting the mean drop error for this subject and visit: %s\n', allSubjects{idx});
        end
    end
end

% estimate correlations between all possible combinations
rho = corr(allMeanDropErrorPerVisit, 'rows', 'complete');
bMask = triu(rho) == 0;
fprintf('Average correlation of mean drop errors across visits: %.3f.\n', mean(rho(bMask)));

% Cronbach's alpha
isValid = all(~isnan(allMeanDropErrorPerVisit), 2);
ca = cronbach(allMeanDropErrorPerVisit(isValid, :)); % cronbach function by Alexandros Leontitsis
fprintf('Cronbach''s alpha: %.3f.\n', ca);
