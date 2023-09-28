%% VIGISTOM TASK 2023 %%
function VIGISTOM(subj_idx)
%% Set path & output structure

global EXP
EXP = [];

EXP.subj_idx = subj_idx;
EXP.nback = 2;

EXP.path2expe = '/Users/lokyan/Library/Mobile Documents/com~apple~CloudDocs/VIGISTOM_TASK/data';

nameOutputFile = strcat(EXP.path2expe,filesep,'Subject_VIGISTOM_structure_',num2str(EXP.subj_idx));
save(nameOutputFile,'EXP');
%% SCREEN

EXP.screenNumber = max(Screen('Screens')); 

EXP.Environment = 'Mac';
EXP.Screen.full = 0;

EXP = NVG_init_screenParameters(EXP);

%% Settings

EXP.Settings.BackgroundColor = [155 155 155];  % Background color in RGB format (from 0 to 255) - greyish
EXP.Settings.FontColor.black    = [0 0 0];        % Letters color in RGB format (from 0 to 255) - black
EXP.Settings.FontColor.red    = [255 0 0]; % red
EXP.Settings.FontColor.white = [255 255 255]; % white
EXP.Settings.FontColor.blue = [0 0 255]; % blue
EXP.Settings.FontColor.green = [0 255 0]; % green
EXP.Settings.FontColor.yellow = [255 255 0]; % yellow

EXP.Settings.LineWidth       = 4;
EXP.Settings.BackgroundColor = EXP.Settings.BackgroundColor./256;
EXP.Settings.FontType       = 'Arial';
EXP.Settings.FontSize       = 40;

%% Initiate Output structure

EXP.LogfileDataMatrix =[];
%% START SCREEN & INITIALIZING PSYCHTOOLBOX
EXP.ifi = 0.0166;
EXP = Screen_start(EXP);

% Wait for trigger from scanner

tStart = KbName('t');

KbQueueCreate();
KbQueueStart();

while true
    [pressed, firstPress] = KbQueueCheck;

    if pressed
        if firstPress(tStart)
            EXP.fMRI_starts = firstPress(tStart);
            break;  % Exit the loop if 't' key is pressed
        end
    end
end

KbQueueStop();
KbQueueRelease();
%% DEFINE TASK PARAMETERS

% Trial duration
EXP.minTrialDuration = 5; % duration of one trial includes the time one series of letters were presented, followed by the presentation of the questions
EXP.maxTrialDuration = 20;

% repetition interval
EXP.minRepInt = 3; % Note that minimum repetition interval cannot be smaller than n-back
EXP.maxRepInt = 6;

EXP.LetterSet= {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'};

% EXP.stim_dur = 1;
EXP.stim_dur = stim_dur;
save(nameOutputFile,'EXP');
%% DEFINE TASK PARAMETERS

EXP.nbackSeries = 0; % counter of number of trials

% initiate KbQueueCheck

KbQueueCreate();
KbQueueStart();

for taskblock = 1:8

    BlockBegin = GetSecs;

    while (GetSecs - BlockBegin) < 90 % Should be 10 minutes (600 seconds), i.e. Run the task until time is up

        % 8 times 10 minutes, put in the number of task block each time

        randLetterSet_ready = 0;

        if randLetterSet_ready == 0

            EXP.TrialDuration = randi([EXP.minTrialDuration, EXP.maxTrialDuration]); % Trial duration here is the length of the letter chain
            EXP.randLetterSet = {};
            EXP = init_jitters(EXP);
            EXP.taskblock = taskblock;
            EXP.nbackSeries = EXP.nbackSeries + 1;

            % log onset of each trial in respect to fMRI start and the
            % duration of each jitter
            EXP.LogfileDataMatrix.onsets.onsetTrial(EXP.nbackSeries) = GetSecs - EXP.fMRI_starts;
            EXP.LogfileDataMatrix.jitters(EXP.nbackSeries, 1:length(EXP.randomIntervals)) = EXP.randomIntervals;

            for i = 1:EXP.TrialDuration
                % Randomly select a letter from LetterSet
                letter = EXP.LetterSet{randi(length(EXP.LetterSet))};
                % Append the selected letter to the string
                EXP.randLetterSet = [EXP.randLetterSet letter];
                msg1 = sprintf('%d-th letter added to randLetterSet string', i);
                disp(msg1);
            end

            for i = 1:(EXP.TrialDuration - nback)

                while strcmp(EXP.randLetterSet{i},EXP.randLetterSet{i+nback})
                    EXP.randLetterSet{i+nback} = EXP.LetterSet{randi(length(EXP.LetterSet))};
                    disp('letter at nback position substituted');
                end

            end

            randLetterSet_ready = 1;

        end

        if randLetterSet_ready == 1

            RepLetterInd = 0; % This variable tracks the index of the letter in EXP.randLetterSet that is currently repeated


            while RepLetterInd <= EXP.TrialDuration

                RepInt = randi([EXP.minRepInt, EXP.maxRepInt]);
                msg2 = sprintf('random Interval: %d', RepInt);
                disp(msg2);

                RepLetterInd = RepInt + RepLetterInd;
                msg3 = sprintf('Index of currently repeated letter: %d', RepLetterInd);
                disp(msg3);

                % if the letter that is to be repeated is within the range of nback
                % + 1 and the length of EXP.randLetterSet:

                if RepLetterInd <= length(EXP.randLetterSet) && (RepLetterInd - nback) > 0

                    EXP.randLetterSet(RepLetterInd) = EXP.randLetterSet(RepLetterInd - nback);
                    disp('nback letter generated');
                    % substitute the 'n' letter before RepLetterInd with the same letter
                    EXP.LogfileDataMatrix.randLetterSet{EXP.nbackSeries, :} = EXP.randLetterSet;

                else
                    % Handle the case when RepLetterInd is out of range
                    disp('EXP.randLetterSet generated');
                    break;
                end

            end

            randLetterSet_ready = 2;
        end

        if randLetterSet_ready == 2

            TargetKey = KbName('1!');
            KeyNext = KbName('4$');

            for iLetter = 1:length(EXP.randLetterSet)

                StimOnsetLogged = 0;

                stim_present = GetSecs;

                button_press = 0;

                while (GetSecs - stim_present) <= (EXP.stim_dur + EXP.randomIntervals(iLetter))

                    if (GetSecs - stim_present) <= EXP.stim_dur

                        letter = EXP.randLetterSet{iLetter};

                        Screen('TextSize', EXP.windowPtr, 80);
                        DrawFormattedText(EXP.windowPtr, letter, 'center', 'center', EXP.Settings.FontColor.white);
                        Screen('Flip', EXP.windowPtr);

                        if StimOnsetLogged == 0

                            StimOnset = GetSecs;
                            StimOnset_fmri = GetSecs - EXP.fMRI_starts;
                            EXP.LogfileDataMatrix.onsets_fmri.letter(EXP.nbackSeries, iLetter) = StimOnset_fmri;

                            StimOnsetLogged = 1;

                        end

                    elseif (GetSecs - stim_present) > EXP.stim_dur

                        Screen('DrawText', EXP.windowPtr, sprintf('+'),EXP.xcenter-13, EXP.ycenter-22,EXP.Settings.FontColor.white);
                        Screen('TextSize',EXP.windowPtr, 80);
                        Screen('Flip', EXP.windowPtr);

                        WaitSecs(EXP.randomIntervals(iLetter));

                    end

                    % acc = 1 -> Hit
                    % acc = -1 -> Miss
                    % acc = 2 -> false alarm
                    % acc = -2 -> correct rejection
                    % acc = 99 -> other keys were pressed

                    if button_press == 0

                        [pressed, firstPress] = KbQueueCheck(); % Check for key presses as long as it is within ITI and as long as nothing is pressed yet

                        if iLetter <= nback

                            if pressed == 0 % correct rejection
                                acc = -2;
                                rt = 0;
                                rt_fmri = 0;

                            elseif pressed && firstPress(TargetKey) % false alarm
                                acc = 2;
                                rt = GetSecs - StimOnset;
                                rt_fmri = GetSecs - EXP.fMRI_starts;
                                button_press = 99;

                            elseif pressed && ~firstPress(TargetKey) % other keys were pressed
                                acc = 99;
                                rt = GetSecs - StimOnset;
                                rt_fmri = GetSecs - EXP.fMRI_starts;
                                button_press = 99;

                            end

                        elseif iLetter > nback

                            if strcmp(EXP.randLetterSet{iLetter}, EXP.randLetterSet{iLetter-nback})

                                if pressed

                                    button_press = 99; % response has already been made, so stop tracking button press

                                    if firstPress(TargetKey) % Hit

                                        acc = 1;
                                        rt = GetSecs - StimOnset;
                                        rt_fmri = GetSecs - EXP.fMRI_starts;

                                    else

                                        acc = 99; % other key was pressed
                                        rt = GetSecs - StimOnset;
                                        rt_fmri = GetSecs - EXP.fMRI_starts;

                                    end

                                else % Miss

                                    acc = -1;
                                    rt = 0;
                                    rt_fmri = 0;

                                end

                            elseif ~strcmp(EXP.randLetterSet{iLetter}, EXP.randLetterSet{iLetter-nback})


                                if pressed == 0 % correct rejection

                                    acc = -2;
                                    rt = 0;
                                    rt_fmri = 0;

                                elseif pressed && firstPress(TargetKey) % false alarm

                                    acc = 2;
                                    rt = GetSecs - StimOnset;
                                    rt_fmri = GetSecs - EXP.fMRI_starts;
                                    button_press = 99;

                                elseif pressed && ~firstPress(TargetKey) % other key was pressed

                                    acc = 99;
                                    rt = GetSecs - StimOnset;
                                    rt_fmri = GetSecs - EXP.fMRI_starts;

                                end

                            end
                        end

                    end

                end

                EXP.LogfileDataMatrix.acc(EXP.nbackSeries,iLetter) = acc;
                EXP.LogfileDataMatrix.rt(EXP.nbackSeries, iLetter) = rt;
                EXP.LogfileDataMatrix.rt_fmri(EXP.nbackSeries, iLetter) = rt_fmri;
                save(nameOutputFile,'EXP');

            end

            EXP = Questionnaire(EXP);

            randLetterSet_ready = 0;

            save(nameOutputFile,'EXP');

        end
    end

    if taskblock ~= 8

    EXP.Break = ['Block ' num2str(taskblock) ' finished. Please take a short break. \n \n'... 
        'When you are ready, press 4 to start the next block.'];

    Screen('TextSize', EXP.windowPtr, 20);
    DrawFormattedText(EXP.windowPtr, EXP.Break, 'Center', 'Center', EXP.Settings.FontColor.white);
    Screen('Flip', EXP.windowPtr); %projects to identified screen

    else

    EXP.End = ['This is the end of the experiment.\n \n'... 
        'Thank you for your participation!'];

    Screen('TextSize', EXP.windowPtr, 20);
    DrawFormattedText(EXP.windowPtr, EXP.End, 'Center', 'Center', EXP.Settings.FontColor.white);
    Screen('Flip', EXP.windowPtr); %projects to identified screen

    end

    response = 0;

    while response == 0
    
        [pressed, firstPress] = KbQueueCheck();

        if pressed

            response = 1;
            
            if firstPress(KeyNext)

                break;
            else
                response = 0;
            end
        end
    end

end

KbQueueStop();
KbQueueRelease();
%% Close Psychtoolbox
sca;