# optimal\_inference\_of\_sameness

This repository is associated with the following paper: 

Van Den Berg, R., Vogel, M., Josić, K., & Ma, W. J. (2012). Optimal inference of sameness. Proceedings of the National Academy of Sciences, 109(8), 3178-3183.

The repository ONLY contains MATLAB code for experiment and data simulation. Please refer to this paper when using the code or data in this repository. If you want model fitting code, please contact Ronald van den Berg through email.

## Data

The data can be downloaded from the Open Science Framework project page, [Exp1](https://osf.io/t4c5m), [Exp2](https://osf.io/mznwy).

The Exp 1 data is one file per set size per subject, while the Exp 2 data is one file per subject, all in .mat files.

Data are saved in a matrix named 'data', where each row is a trial and the columns are as follows:

- Column 1: mean of the distribution over stimulus orientations for this trial
- Column 2: sample standard deviation of the distribution over stimulus orientations (0 in "same" trials, varies in "different" trials)
- Column 3: population standard deviation of the distribution over stimulus orientations (0 in "same" trials, 10 in "different" trials)
- Column 4: response (0="same", 1="different")
- Column 5: response correctness (0=incorrect, 1=correct)
- Column 6: response time
- Column 7: stimulus time (as measured)
- Columns 8-15: reliability values for the stimuli (only first N are used, where N is the set size on this trial)
- Columns 16-23: stimulus orientations (only first N are used, where N is the set size on this trial)
- Columns 24-31: stimulus locations (measured as the angle of the polar coordinate representation; only first N are used, where N is the set size on this trial)

