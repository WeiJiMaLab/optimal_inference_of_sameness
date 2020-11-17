# optimal\_inference\_of\_sameness

This repository is associated with the following paper: 

Van Den Berg, R., Vogel, M., Josić, K., & Ma, W. J. (2012). Optimal inference of sameness. Proceedings of the National Academy of Sciences, 109(8), 3178-3183.

The repository only contain MATLAB code for experiment and data simulation (to-update). Please refer to this paper when using the code in this repository.

## Data (to-update)

The data can be downloaded from the Open Science Framework project page (link-to-update).

The Exp 1 data is one file per set size per subject, while the Exp 2 data is one file per dubject, all in .mat files.

Data are saved in a variable named 'data' that contains trial-by-trial informaion as

data(trial_number,:) = [setMu setSigma realSetSigma RESPIDX CORRECT RT REALSTIMTIME setEpsilons setOrts setPosThetas];

where each column is

1. setMu = mean of the distribution over stimulus orientations for this trial
2. setSigma = sample standard deviation of the distribution over stimulus orientations (0 in "same" trials, varies in "different" trials)
3. realSetSigma = population standard deviation of the distribution over stimulus orientations (0 in "same" trials, 10 in "different" trials if i remember correctly)
4. respIDX = response (should be 2 possible values, where one of them is "same" the other "different")
5. CORRECT = probably 0 when response was incorrect and 1 otherwise
6. RT = response time
7. REALSTIMTIME = stimulus presentation time (just for check)
8. setEpsilons = 8 columns with reliability values for the stimuli (only first N are used, where N is the set size on this trial)
9. setOrts = 8 columns with stimulus orientations (only first N are used, where N is the set size on this trial)
10. setPosThetas = 8 columns with stimulus locations (measured as the angle of the polar coordinate representation of this position; only first N are used, where N is the set size on this trial)

