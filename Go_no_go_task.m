%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programming Project % Go-nogo Task%
% WiSe 22/23 %%%%%%%%%%%%%%%%%%%%%%%%
% Created by: %%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                                                                                                                          
% Zeynep, Janice, Lokyan & Nele %%%%%                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                                                                                                                                                                            

%Creating the function
function Go_no_go_task

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HOUSEKEEPING! HOUSEKEEPING!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

%Clear Matlab/Octave window:
clc;
close all; clearvars;

% create paths to our folders + constants to add to subject ID numbers to
% avoid overlapping

% dir_zep = "C:\Users\Zeynep Elif Sab\Desktop\GO_NOGO_EXP"; 
%zeynep = 100;

% dir_janice = '/Users/apple/Documents/MATLAB/data'; 
%janice = 200;

dir_lokyan = '/Users/lokyan/Documents/Task_Programming/'; 
lokyan = 300;

% dir_nele = '/Users/nelebehrens/Desktop/no_go:data:';
% nele = 400;

% Insert your name here to go to your own directory
cd(dir_nele);
              
%% participant ID
resultsfolder = './data/' % go to the folder 'data' in the directory
numofsubj = length(dir(resultsfolder)) - 2 
subj = numofsubj + lokyan % insert your name here

%% Open a text file for each participant to store results

% Set column names
colHeaders = ["subj", "block", "trial number", "phase label", "stimulus type", "accuracy","reaction time"];

% open a text file for the participant under the resultsfolder
datafilename = strcat(resultsfolder,'GoNogo_',num2str(subj),'.txt'); 

% write column names into the text file
fopen(datafilename,'wt');


%% Prepare variables and materials for the experiment

% define number of trials
ntrials = 150;
% define stimulus presentation time
stim_dur = 0.5;
% define fixation cross presentation time
fixcross_dur = 0.8;

% counterbalance condition_list: Two laptops would start with 10 (happy-go) block
% while the other two laptops would start with 20 (sad-go) block
% condition_list = [10 20]; % JANICE & ZEYNEP
condition_list = [20 10]; % LOKYAN & NELE

% define our results matrix size - results would be stored into this matrix
% and written into the text file by the end of every trial
results = NaN * ones(length(ntrials*2),length(colHeaders));

% prepare stimulus list
for block = 1:2 % number of blocks 
    stimuli_list(block) = {cell(1,ntrials)}; % each block should contain a 1x150 cell array
    if condition_list(block) == 10
        emotion_presented = [repmat("happy",1,ntrials*0.8) repmat("sad",1,ntrials*0.2)]; %add 80% go (happy) and 20% no-go (sad) stimuli
    elseif condition_list(block) == 20
        emotion_presented = [repmat("sad",1,ntrials*0.8) repmat("happy",1,ntrials*0.2)]; %add 80% go (sad) and 20% no-go (happy) stimuli
    end  

    stimuli_list{block} =  emotion_presented(randperm(ntrials)); %shuffle emotions that are presented in one block
end

disp('Stimulus list has been prepared');

% prepare instructions
Intro_1 = ['In this experiment you will be presented with two emoji faces:\n \n' ...
              ' (1) A smiling emoji and (2) a frowning emoji. \n \n' ...
              '  You will be asked to press the spacebar AS QUICKLY AS POSSIBLE as soon as an emoji face appears. \n \n' ...
              '  (Press the space bar to continue)' ];

Intro_2 = ['You will now be asked to push the spacebar every time the smiling emoji appears. \n \n' ...
    'DO NOT press the spacebar when you see the frowning emoji! \n \n'...
             '  (Press the spacebar to begin)' ];

Intro_3 = ['You will now be asked to push the spacebar every time the frowning emoji appears. \n \n' ...
    'DO NOT press the spacebar when you see the smiling emoji! \n \n'...
             '  (Press the spacebar to begin)' ];

% prepare fixation cross
    FixCr=ones(20,20)*255;
    FixCr(10:11,:)=0;
    FixCr(:,10:11)=0;

    %% Launch PTB %%%
    try
        %Optional, run if your laptop has problems with PTB
        Screen('Preference','SkipSyncTests',1);

        AssertOpenGL;
        KbName('UnifyKeyNames');

        %set up screen
        Screens = Screen('Screens'); % identify screens
        screenNumber = max(Screens);

        [myscreen,rect] = Screen('OpenWindow',screenNumber); %% DELETE 3RD AND 4TH ARGUMENT WHEN THE SCRIPT IS FINALISED

        %hide cursor for aesthetic reasons
        HideCursor;

        % set target key for participant's response
        TargetKey = KbName('space');

        % block_order = ["happy" "sad"]; % for JANICE & ZEYNEP
        block_order = ["sad" "happy"]; % for LOKYAN & NELE

        fixcross = Screen('MakeTexture',myscreen,FixCr); % save the fixation cross as a texture (for displying later)

        while KbCheck; end

        %1. Introduction
        Screen('TextSize', myscreen, 16);
        DrawFormattedText(myscreen, Intro_1, 'center', 'center');
        Screen('Flip', myscreen); %projects to identified screen
        KbWait([], 3); %waits for keyboard press before continuing

        %add second intro to explain which the go-stimuli is --> switch to Intro_3, if we start with sad as go
        Screen('TextSize', myscreen, 16);
        DrawFormattedText(myscreen, Intro_3, 'center', 'center');
        Screen('Flip', myscreen); %projects to identified screen
        KbWait([], 3); %waits for keyboard press before continuing

        rt = zeros(2,ntrials); % open a list for rt
        % By wrapping cell(1,ntrials) in curly brackets, we are creating nested
        % cells inside the cell for 'block'
        acc = zeros(2,ntrials); % open a list for acc


        %% Presentation of stimuli

        % prepare row_number counter so we can go over 10 trials in the results (b/c 2 blocks * 10 trials)
        row_number = 0;

        for block = 1:2

            % go_stim is the stimulus to which the participants should respond to i.e. happy or sad
            % This variable would be used in the loop below
            go_stim = block_order(block);

            if block == 2 %add intro between blocks --> switch to intro_2 for happy as go in second block
                Screen('TextSize', myscreen, 16);
                DrawFormattedText(myscreen, Intro_2, 'center', 'center');
                Screen('Flip', myscreen); %projects to identified screen
                KbWait([], 3); %waits for keyboard press before continuing
            end

            Screen('DrawTexture', myscreen, fixcross); % draw fixation cross to backbuffer
            Screen('Flip', myscreen) %projects to identified screen
            WaitSecs(fixcross_dur); %wait until continue

            for trial = 1:length(stimuli_list{1,block}) % LENGTH = TOTAL NUMBER OF ENTRIES IN THE LIST FOR THE BLOCK OF INTEREST
                
                row_number = row_number + 1;
                
                randomjitter = randi([-80 80], 1)/1000;
                ITI = (stim_dur + fixcross_dur) + randomjitter; 
                % ITI is the time period from the start of one face to the start of the next face

                msg = sprintf('%s - Go Block. Trial identity: %s',go_stim,stimuli_list{1,block}(trial));
                disp(msg);

                button_press = 0; 
                % button_press allows us to stop KbCheck as soon as a key is pressed.
                % If nothing is pressed, KbCheck would keep running until time is
                % up (ITI is reached).
                stimonset_register = 0; 
                % stimonset_register allows stimulus onset to be registered
                % ONE time only (i.e., the first time the stimulus is presented), 
                % even though the stimulus would be presented over
                % and over again, as long as it is within the stimulus
                % presentation period (stim_dur)

                trial_start = GetSecs; % The trial starts now, time onset is saved into trial_start

                while (GetSecs - trial_start) <= ITI % if current time - trial onset <= ITI, we stay in the while loop

                    if (GetSecs - trial_start) <= stim_dur % if current time - trial onset <= stim_dur, present face stimulus

                        filename = strcat(stimuli_list{1,block}(trial), 'face.jpeg'); % stimuli_list{1,block}(trial) RETURNS "HAPPY" OR "SAD"

                        % This is to allow the info about the stimulus type to be written into the results matrix of double type
                        if strcmp(stimuli_list{1,block}(trial),"happy")                            
                            stimulustype = 88; % 88 is the code for happy
                        else
                            stimulustype = 44; % 44 is the code for sad
                        end

                        stim = imread(filename); % read the .jpeg file and store it in a matrix
                        tex=Screen('MakeTexture', myscreen, stim); % save the stimulus as a texture
                        Screen('DrawTexture', myscreen, tex); % display the stimulus
                        Screen('Flip', myscreen);
                        if stimonset_register == 0 
                            stimonset = GetSecs; % register stimulus onset
                            stimonset_register = 99; % so that stimonset can only be registered once
                        end

                    elseif (GetSecs - trial_start) > stim_dur % if current time - trial onset > stim_dur, present fixation cross

                        Screen('DrawTexture', myscreen, fixcross);% prepare fixation cross
                        Screen('Flip', myscreen); % display the fixation cross

                    end

                    [KeyIsDown,secs,keyCode] = KbCheck; % Check for key presses as long as it is within ITI

                    if button_press == 0 

                        if KeyIsDown == 1

                            disp('some key was pressed');
                            rt(block, trial) = (secs-stimonset)*1000; % SAVE REACTION TIME
                            RT = (secs-stimonset)*1000;
                            msg = sprintf('response time: %d ,stimulus onset: %d ,reaction time: %d', secs, stimonset, RT);
                            disp(msg);

                            button_press = 99; % KbCheck would not be perforced anymore in this trial.

                            if strcmp(stimuli_list{1,block}(trial), go_stim) % IF TRIAL IN STIMULUS LIST (E.G. HAPPY) = GO-TRIAL (E.G HAPPY)

                                if keyCode(TargetKey) %%% IF THE KEY PRESS IS THE SPACEBAR

                                    disp('Spacebar was down');
                                    acc(block, trial) = 1; % SAVE ACCURACY = 1
                                    disp('Accuracy = 1 was saved in the logfile');                                   

                                else %%% IF THE KEY PRESS IS SOMETHING OTHER THAN THE SPACEBAR

                                    disp('Something other than spacebar was down');
                                    acc(block, trial) = 0; % SAVE ACCURARY = 0
                                    disp('Accuracy = 0 was saved in the logfile');                                   

                                end

                            else % IF TRIAL IN STIMULUS LIST (E.G. HAPPY) = NOGO-TRIAL (E.G SAD)

                                disp('something was pressed');
                                acc(block, trial) = 0; % SAVE ACCURACY = 0
                                disp('Accuracy = 0 was saved in the logfile');        

                            end

                        else % If no key was pressed

                            button_press = 0; % let KbCheck check for key presses
                                % as much as possible within ITI for any potential
                                % presses

                            if strcmp(stimuli_list{1,block}(trial), go_stim) % IF TRIAL IN STIMULUS LIST (E.G. HAPPY) = GO-TRIAL (E.G HAPPY)

                                disp('nothing was pressed');
                                acc(block, trial) = 0; % SAVE ACCURACY = 0
                                disp('Accuracy = 0 was saved in the logfile');                            

                            else % IF TRIAL IN STIMULUS LIST (E.G. HAPPY) = NOGO-TRIAL (E.G HAPPY)

                                disp('nothing was pressed');
                                acc(block, trial) = 1; % SAVE ACCURARY = 1
                                disp('Accuracy = 1 was saved in the logfile');

                            end

                            %Accuracy = 1 --> correct, Accuracy = 0 -->incorrect

                        end
                    end
                end

                % enter results in matrix
                results(row_number,:) = [subj, block, trial, stimulustype, condition_list(block), acc(block,trial), rt(block,trial)] %[phaselabel];%, stimuli_list{trial}, trialno, acc, rt];

                % write results to comma delimited text file 
                dlmwrite(datafilename, results, 'delimiter', ',', 'precision', 6);
            
                % save results matrix to a .mat file to load as a variable
                % during data analysis
                save(strcat('var_GoNogo',num2str(subj),'.mat'),'results')
            end
        end

        %calculate and display performance feedback for all
        accuracy = 0;
        for i = 1:size(acc,1) % per row
            for j = 1:size(acc,2) % per column
                if acc(i,j) == 1
                    accuracy = accuracy + 1;
                end
            end
        end

        correct_perc = accuracy/(ntrials*2)

        %calculate and display performance feedback for separate blocks
        accuracy_block1 = 0
        accuracy_block2 = 0
        row1 = acc(1,:)
        row2 = acc(2,:)

        for i = 1:size(row1,1) % per row
            for j = 1:size(row1,2) % per column
                if row1(i,j) == 1
                    accuracy_block1 = accuracy_block1 + 1;
                end
            end
        end

        for i = 1:size(row2,1) % per row
            for j = 1:size(row2,2) % per column
                if row2(i,j) == 1
                    accuracy_block2 = accuracy_block2 + 1;
                end
            end
        end

        correct_perc_block1 = accuracy_block1/ntrials;
        correct_perc_block2 = accuracy_block2/ntrials;

        % write the outro message
        out_message = ['The experiment has ended. \n \n' ...
            'You were correct in ' num2str(correct_perc_block1*100) '% of trials in the first block, \n \n'...
            'while you were correct in ' num2str(correct_perc_block2*100) '% of trials in the second block.\n \n'...
            'Thank you for your compulsory participation. \n \n' ...
            'NOW PRESS SPACEBAR TO EXIT.' ];

        % load the buffer
        DrawFormattedText(myscreen, out_message, 'center', 'center');

        % Update the display to show the instruction text:
        Screen('Flip', myscreen);

        %wait for keystroke
        KbWait([], 3);

        sca;
        ShowCursor;

        Screen('CloseAll');
    catch
        sca;
        ShowCursor;
        psychrethrow(psychlasterror)
    end