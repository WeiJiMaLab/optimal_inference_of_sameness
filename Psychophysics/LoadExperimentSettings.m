function settings = LoadExperimentSettings(subjid, condNr)

settings = [];

% condition-independent settings
settings.pSame      =         .5;  % prior probability of having a "same" trial
settings.ecc        =          8;  % eccentricity at which stimuli are shown (in deg visual angle)
settings.ITT        =        500;  % inter-trial-time (msec)
settings.bglum      =         30;  % background luminance (cd/m^2)
settings.fglum      =         40;  % foreground luminance (cd/m^2)
settings.stimtime   =         55;  % stimulus time (msec)
settings.sigmaNZ    =         10;
settings.stimsize   =         .5;  % stimulus size, in deg^2
settings.posJitterSigma =    .25;
settings.noTrials    =       150;  % number of trials 

settings.bglum      =         10;  % background luminance (cd/m^2)
settings.fglum      =         50;  % foreground luminance (cd/m^2)


if condNr==0

    % CONDITION 0: PRACTICE
    
    settings.exptype     = 'practice';  % training condition, meaning that stimulus time will depend on performance
    settings.epsilons    =         .94;
    settings.multiEps   =           0;  % if 0, all stimuli within a trial have same epsilon (drawn from settings.epsilons)
    settings.N           =          6;  % vector with all possible numbers of stimuli
    settings.stimtime    =        300; 
    settings.minstimtime =         50;
    settings.maxstimtime =       1000;
    settings.feedback    =          1;
    settings.nPractice   =          0;
    
elseif condNr==1

    % CONDITION 1: SIGMA_NZ THRESHOLDS
    
    settings.exptype     =   'sigma_threshold';  % threshold measurement
    settings.sigmaNZ     =   linspace(5,15,12);  % set std. dev. in "different" trials (if this is a vector, a random sigma_NZ will be drawn from this, on "different" trials)
    settings.epsilons =            .8;
    settings.multiEps   =           0;  % if 0, all stimuli within a trial have same epsilon (drawn from settings.epsilons)
    settings.N           =        [4];  % vector with all possible numbers of stimuli
    settings.feedback    =          1;  % feedback yes/no
    settings.nPractice   =          0;  % #practice trials preceding the 'actual' trials
    
elseif condNr==2
    
    % CONDITION 2: ELLIPSE ECCENTRICITY THRESHOLDS

    settings.exptype     =   'eps_threshold';  % threshold measurement
    settings.epsilons    =   linspace(.5,.94,10);
    settings.multiEps    =          0;  % if 0, all stimuli within a trial have same epsilon (drawn from settings.epsilons)
    settings.N           =          6;  % vector with all possible numbers of stimuli
    settings.feedback    =          1;  % feedback yes/no
    settings.nPractice   =          0;  % #practice trials preceding the 'actual' trials
    
elseif condNr==3 
    
    % CONDITION 3: ACTUAL EXPERIMENT, LOW EPSILON    

    settings.exptype    = 'actual_low';  % threshold measurement
    settings.multiEps  =           0;  % if 0, all stimuli within a trial have same epsilon (drawn from settings.epsilons)
    settings.N          =          6;  % vector with all possible numbers of stimuli
    settings.feedback   =          1;  % feedback yes/no
    settings.nPractice  =          0;  % #practice trials preceding the 'actual' trials

    % compute epsilon_low based on results from condition 2:
    [eps_low eps_high] = compute_epsilon_thresholds(subjid,0);
    settings.epsilons = eps_low;
    
elseif condNr==4 
    
    % CONDITION 4: ACTUAL EXPERIMENT, MIXED EPSILON

    settings.exptype    = 'actual_mixed';  % threshold measurement
    settings.multiEps   =          1;  % if 0, all stimuli within a trial have same epsilon (drawn from settings.epsilons)
    settings.N          =          6;  % vector with all possible numbers of stimuli
    settings.feedback   =          1;  % feedback yes/no
    settings.nPractice  =          0;  % #practice trials preceding the 'actual' trials

    % compute epsilon_low based on results from condition 2:
    [eps_low eps_high] = compute_epsilon_thresholds(subjid,0);
    settings.epsilons = [eps_low eps_high];
    
elseif condNr==5 
    
    % CONDITION 5: ACTUAL EXPERIMENT, HIGH EPSILON

    settings.exptype    = 'actual_high';  % threshold measurement
    settings.multiEps  =           0;  % if 0, all stimuli within a trial have same epsilon (drawn from settings.epsilons)
    settings.N          =          6;  % vector with all possible numbers of stimuli
    settings.feedback   =          1;  % feedback yes/no
    settings.nPractice  =          0;  % #practice trials preceding the 'actual' trials

    % compute epsilon_low based on results from condition 2:
    settings.epsilons = [.94];
   
%-%- EXPERIMENT 2 %-%-%    

elseif condNr==6
    
    % CONDITION 6: ACTUAL EXPERIMENT, HIGH EPSILON, N=6

    settings.exptype    =  'exp2_n8';  % threshold measurement
    settings.multiEps   =          0;  % if 0, all stimuli within a trial have same epsilon (drawn from settings.epsilons)
    settings.N          =          8;  % vector with all possible numbers of stimuli
    settings.feedback   =          1;  % feedback yes/no
    settings.nPractice  =          0;  % #practice trials preceding the 'actual' trials
    settings.epsilons   =      [.94];

elseif condNr==7
    
    % CONDITION 7: ACTUAL EXPERIMENT, HIGH EPSILON, N=4

    settings.exptype    =  'exp2_n4';  % threshold measurement
    settings.multiEps  =           0;  % if 0, all stimuli within a trial have same epsilon (drawn from settings.epsilons)
    settings.N          =          4;  % vector with all possible numbers of stimuli
    settings.feedback   =          1;  % feedback yes/no
    settings.nPractice  =          0;  % #practice trials preceding the 'actual' trials
    settings.epsilons   =      [.94];

elseif condNr==8
    
    % CONDITION 8: ACTUAL EXPERIMENT, HIGH EPSILON, N=2

    settings.exptype    =  'exp2_n2';  % threshold measurement
    settings.multiEps  =           0;  % if 0, all stimuli within a trial have same epsilon (drawn from settings.epsilons)
    settings.N          =          2;  % vector with all possible numbers of stimuli
    settings.feedback   =          1;  % feedback yes/no
    settings.nPractice  =          0;  % #practice trials preceding the 'actual' trials
    settings.epsilons   =      [.94];
    
elseif condNr==9
    
    % CONDITION 9: COLOR EXPERIMENT 
    settings.sigmaNZ    =  10;
    settings.exptype    =  'exp3_color';  % color experiment
    settings.N          =  [2 4 8];       % vector with all possible numbers of stimuli
    settings.feedback   =  1;             % feedback yes/no
    settings.epsilons   =  1;             % doesnt have a meaning in this experiment
    settings.multiEps   =  0;             % doesnt have a meaning in this experiment
    settings.ecc        =  6;  % eccentricity at which stimuli are shown (in deg visual angle)
    settings.stimsize   =  .4;  % stimulus size, in deg^2
    settings.posJitterSigma =  0;

elseif condNr==10
    
    % CONDITION 10: COLOR EXPERIMENT PRACTICE
    settings.sigmaNZ    =  10;
    settings.exptype    =  'exp3_color';  % color experiment
    settings.N          =  [2 4 8];       % vector with all possible numbers of stimuli
    settings.feedback   =  1;             % feedback yes/no
    settings.epsilons   =  1;             % doesnt have a meaning in this experiment
    settings.multiEps   =  0;             % doesnt have a meaning in this experiment
    settings.ecc        =  6;  % eccentricity at which stimuli are shown (in deg visual angle)
    settings.stimsize   =  .4;  % stimulus size, in deg^2
    settings.posJitterSigma =  0;
    settings.noTrials   = 25;
elseif condNr==11
    
    % CONDITION 11: DEMO
    settings.stimtime   =  500;  % stimulus time (msec)
    settings.sigmaNZ    =  10;
    settings.exptype    =  'exp3_color';  % color experiment
    settings.N          =  [2 4 8];       % vector with all possible numbers of stimuli
    settings.feedback   =  1;             % feedback yes/no
    settings.epsilons   =  1;             % doesnt have a meaning in this experiment
    settings.multiEps   =  0;             % doesnt have a meaning in this experiment
    settings.ecc        =  6;  % eccentricity at which stimuli are shown (in deg visual angle)
    settings.stimsize   =  .4;  % stimulus size, in deg^2
    settings.posJitterSigma =  0;
    settings.noTrials   = 6;
else
       
    fprintf('FAIL!');
    
end


if (settings.sigmaNZ == -1)
    fprintf('\n\n**************************************************\n');
    fprintf('** Could not determine sigmaNZ for this subject **\n');
    fprintf('**************************************************\n\n');
    error('TERMINATED');
end

if ~isempty(intersect(settings.epsilons,-1))
    fprintf('\n\n************************************************************\n');
    fprintf(    '** Could not determine epsilon threshold for this subject **\n');
    fprintf(    '************************************************************\n\n');
    error('TERMINATED');
end
