%==========================================================================
% This script preprocesses the behavioral data.
%
% Lukas Kunz, 2023
%==========================================================================

% start
clear; clc; close all;

% paths
paths       = [];
paths.beh   = 'G:\My Drive\Manuscripts\2023_Boecker_SportsAndSpatialMemory\Behavior\';
paths.save  = 'G:\My Drive\Manuscripts\2023_Boecker_SportsAndSpatialMemory\Preprocessing_20230319\';

% settings
param                   = [];
param.direc.maxAbsYaw   = 32768;

% logfiles
logFiles    = dir(strcat(paths.beh, '*.log'));

%% loop through subjects
for iSub = 1:size(logFiles, 1)
    
    % subject name
    subjName        = strrep(logFiles(iSub).name, '.log', '');

    % report
    fprintf('\n\n==== SUBJECT SESSION: %s.\n', subjName);
    
    % save path for this subject
    thisSavePath    = strcat(paths.save, subjName, filesep);
    mkdir(thisSavePath);
    
    % skip this subject, if the final output exists already
    if exist(strcat(thisSavePath, 'behInfo.mat'), 'file') > 0
        continue;
    end

    %% original behavioral data
    
    % file handle
    fid             = fopen(fullfile(logFiles(iSub).folder, logFiles(iSub).name));
    textdata        = [];
    
    % read logfile line-by-line
    iline = 1;
    tline = fgetl(fid);
    while ischar(tline)
        
        % combine all text lines into one variable
        textdata{iline, 1} = tline;
        
        % increase line index
        iline = iline + 1;
        
        % read next line
        tline = fgetl(fid);
    end
    
    % close logfile
    fclose(fid);
    
    % report progress
    fprintf('You have read %d text lines...\n', iline);

    % extract relevant data
    origData    = textdata(contains(textdata, 'ScriptLog'), :);
    
    %% ==================================================================== process trial information ==> trialInfo
    
    % number of trials the subject completed
    numTrials   = sum(contains(origData, 'TrialPhase2Counter ='));
    
    % report
    fprintf('\nExtracting relevant trial information ==> "trialInfo.mat".\n');
    fprintf('Number of trials the subject completed: %d.\n', numTrials);
    
    % create "trials" variable with information for: trialidx, object,
    % ITI-onset, cue-onset, retrieval-onset, feedback-onset,
    % reencoding-onset, grab-onset, x-correct, y-correct, x-drop, y-drop,
    % drop error
    trials      = nan(numTrials, 13);
    bTestPhase  = false; % whether the test phase has begun
    for iLine = 1:size(origData, 1)
        
        % information from this line
        thisLineText        = origData{iLine, :};
        thisLineSplits      = strsplit(thisLineText, ' ');
        if iLine < size(origData, 1)
            nextLineText    = origData{iLine + 1, :};
            nextLineSplits  = strsplit(nextLineText, ' ');
        end
        
        % extract information
        if regexp(thisLineText, 'TrialPhase2Counter =')
            
            % trial index
            trialIdx                = str2double(thisLineSplits{4});
            trials(trialIdx, 1)     = trialIdx;
            
            % note that the test phase has begun
            bTestPhase              = true;
            
        elseif regexp(thisLineText, ' Cue ')
            
            % object identity
            trials(trialIdx, 2)     = str2double(thisLineSplits{5});
        
        elseif regexp(thisLineText, 'ITI_start')
            
            % timepoint of ITI start
            trials(trialIdx, 3)     = str2double(thisLineSplits{3});

        elseif regexp(thisLineText, 'CUE_start')
            
            % timepoint of cue start
            trials(trialIdx, 4)     = str2double(thisLineSplits{3});
            
        elseif regexp(thisLineText, 'MyExp2Spawner.bReadyDrop True')
            
            % timepoint of retrieval start
            trials(trialIdx, 5)     = str2double(nextLineSplits{3}); % cave: timepoint is in the next line
            
        elseif regexp(thisLineText, ' Drop ')
            
            % response location (x/y)
            trials(trialIdx, [11, 12])  = [str2double(thisLineSplits{7}), str2double(thisLineSplits{9})];
            
        elseif regexp(thisLineText, 'DropBeep')
            
            % timepoint of feedback
            trials(trialIdx, 6)     = str2double(thisLineSplits{2});
            
        elseif regexp(thisLineText, ' Show ')
            
            % skip if test phase has not begun yet
            if bTestPhase == false
                continue;
            end
            
            % correct location (x/y)
            trials(trialIdx, [9, 10])   = [str2double(thisLineSplits{7}), str2double(thisLineSplits{9})];
            
        elseif regexp(thisLineText, 'how accurately placed')
            
            % drop error (computed by Unreal)
            trials(trialIdx, 13)    = str2double(thisLineSplits{6});
            
        elseif regexp(thisLineText, 'STOP_SMILEYFB')
            
            % timepoint of re-encoding start
            trials(trialIdx, 7)     = str2double(thisLineSplits{4});
            
        elseif regexp(thisLineText, 'GrabBeep')
            
            % timepoint of grab
            trials(trialIdx, 8)     = str2double(thisLineSplits{2});
        end
    end
    
    % sanity check that time points are strictly increasing
    tmpTimes    = transpose(trials(:, 3:8));
    if any(diff(tmpTimes(:)) < 0)
        error('Timepoints in "trials" are not strictly increasing.');
    end

    %% convert into "trialInfo" table and save
    
    % variable names for table
    varNames    = {'TrialIdx', 'Object', 'ITI', 'Cue', 'Retrieval', 'Feedback', 'Reencoding', 'Grab', 'xCorrect', 'yCorrect', 'xResponse', 'yResponse', 'DropError'};
    
    % convert data into table
    trialInfo   = array2table(trials, 'VariableNames', varNames);
        
    % save trial information
    save(strcat(thisSavePath, 'trialInfo'), 'trialInfo');
    
    %% ==================================================================== process behavioral information ==> behInfo
    
    % report
    fprintf('\nExtracting relevant behavioral information ==> "behInfo.mat".\n');
    
    % preallocate output for relevant data
    behData  	= nan(size(origData, 1), 7); % time, x, y, z, yawUnreal, trial idx, trial phase
    
    % initialize trial index and trial phase
    trialIdx    = nan;
    trialPhase  = nan;
    
    % loop through text lines
    for iLine = 1:size(origData, 1) - 1
        
        % information from this and the next line
        thisLineText = origData{iLine, :};
        nextLineText = origData{iLine + 1, :}; % needed for yawUnreal
        
        % update current trial index
        if regexp(thisLineText, ' TrialPhase2Counter ')
            splits    	= strsplit(thisLineText, ' ');
            trialIdx   	= str2double(splits{4});
        end
        
        % update current trial phase
        if regexp(thisLineText, 'ITI_start')
            trialPhase 	= 1; % ITI
        elseif regexp(thisLineText, 'CUE_start')
            trialPhase 	= 2; % cue
        elseif regexp(thisLineText, 'bReadyDrop True')
            trialPhase 	= 3; % retrieval
        elseif regexp(thisLineText, 'DropBeep')
            trialPhase 	= 4; % feedback
        elseif regexp(thisLineText, 'STOP_SMILEYFB')
            trialPhase 	= 5; % re-encoding
        elseif regexp(thisLineText, 'GrabBeep')
            trialPhase 	= 6; % grab
        end
        
        % get additional information and collect all information whenever
        % subject is at a specific location
        if regexp(thisLineText, 'Location X')
            
            % time and xyz-location information
            splits              = strsplit(thisLineText, ' ');
            behData(iLine, 1:4) = [str2double(splits{3}), str2double(splits{6}), str2double(splits{8}), str2double(splits{10})];
            
            % heading information (cave: contained in the next line)
            splits              = strsplit(nextLineText, ' ');
            behData(iLine, 5)   = str2double(splits{10});
            if ~contains(nextLineText, 'Yaw')
                error('Next line does not contain "Yaw" information.'); % sanity check
            end
            
            % trial-index information
            behData(iLine, 6)   = trialIdx;
            
            % trial-phase information
            behData(iLine, 7)   = trialPhase;
        end
    end
    
    % remove lines that do not have a time stamp
    fprintf('Size of "behData" before cutting lines without a time stamp: \t%d x %d.\n', size(behData));
    behData = behData(~isnan(behData(:, 1)), :);
    fprintf('Size of "behData" after cutting lines without a time stamp: \t%d x %d.\n', size(behData));
    
    %% convert into "behInfo" table and add additional variables
    
    % variable names for the table
    varNames    = {'time', 'x', 'y', 'z', 'yawUnreal', 'trialIdx', 'trialPhase'};
    
    % convert data into table
    behInfo     = array2table(behData, 'VariableNames', varNames);
    
    % durations
    behInfo.durations       = [diff(behInfo.time); median(diff(behInfo.time), 'omitnan')]; % add a median duration at the end

    % translational distances and speed
    behInfo.distances       = [sqrt(diff(behInfo.x) .^ 2 + diff(behInfo.y) .^ 2); 0];
    behInfo.speed  	        = behInfo.distances ./ behInfo.durations;
    
    % angular distances and speed
    behInfo.yaw             = (behInfo.yawUnreal ./ param.direc.maxAbsYaw) .* pi; % yaws in radians
    behInfo.yawDistances    = abs([angdiff(behInfo.yaw); 0]);
    behInfo.yawSpeed        = behInfo.yawDistances ./ behInfo.durations;
    
    %% save output
    
    % save behavioral information
    save(strcat(thisSavePath, '\behInfo'), 'behInfo');
end

% report
fprintf('\nBehavioral processing completed\n.');