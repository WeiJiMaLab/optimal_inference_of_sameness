% function data = generate_fake_data(sigma_low,sigma_high,psame,sigma_s,N,nTrials)
% 
% Generate synthetic data from the Bayesian observer model for the experiments
% described in "Optimal inference of sameness" by Van den Berg, Vogel, Josic, 
% and Ma, PNAS 2012.
% 
% INPUT
%  sigma_low   : standard deviation of the noise distribution corresponding
%                to the stimuli with low reliability
%  sigma_high  : standard deviation of the noise distribution corresponding
%                to the stimuli with high reliability
%  psame       : subject prior for "same" displays (note that this only affects
%                the inference; the prior in the generative model is fixed to 0.5)
%  sigma_s     : subject's estimate of sigma_s (again, not that in the
%                generative model this value is fixed to 10, as in the experiment)
%  N           : set size 
%  nTrials     : number of synthetic trials to generate
%
% OUTPUT
%  data.mu     : means of the distributions from which the sets of stimuli were drawn
%  data.sigma  : std of the distributions from which the sets of stimuli were drawn
%  data.C      : sameness values (-1=different trial, 1=same trial)
%  data.C_hat  : responses (-1=different trial, 1=same trial)
%  data.stimulus_matrix: nTrials x N matrix with stimulus values 
%  data.reliability_matrix: nTrials x N matrix with reliability values
%
% Example:
%   data = generate_fake_data(3,4,0.5,10,4,500);

% Written by RvdB, Feb 2012

function data = generate_fake_data(sigma_low,sigma_high,psame,sigma_s,N,nTrials)

% draw stimuli
data.mu = rand(nTrials,1)*180;         % draw random mean orientation for each trial
data.C = sign(rand(nTrials,1)-0.5);    % C=-1: different trial; C=1: same trial;
data.sigma = (data.C==-1)*10 + (data.C==1)*0; % set sigma of stimulus distribution for each trial (0 on "same" trials, 10 on "different" trials)
data.stimulus_matrix = normrnd(repmat(data.mu,1,N),repmat(data.sigma,1,N)); % draw stimuli for all trials at once
data.reliability_matrix = rand(nTrials,N)>.5; % set stimuli randomly to low/high reliability

% draw internal representations
sigma_int_mat = (data.reliability_matrix==0)*sigma_low + (data.reliability_matrix==1)*sigma_high;
x = normrnd(data.stimulus_matrix,sigma_int_mat);

% compute d and generate responses
w = 1./sigma_int_mat.^2;
w_tilde = 1./(sigma_int_mat.^2+sigma_s.^2);
for ii=1:nTrials
    A = w_tilde(ii,:)'*w_tilde(ii,:)/sum(w_tilde(ii,:)) - w(ii,:)'*w(ii,:)/sum(w(ii,:)) + diag(w(ii,:)-w_tilde(ii,:));
    d(ii) = 0.5 * (-x(ii,:)*A*x(ii,:)' + sum(log(w(ii,:)./w_tilde(ii,:))) - log(sum(w(ii,:))/sum(w_tilde(ii,:)))) + log(psame/(1-psame));
end
data.C_hat = sign(d');  % C_hat=-1: "different" response; C_hat=+1: "same" response

fprintf('Generated a dataset with %d synthetic trials.\nPerformance is %2.1f%% correct\n',nTrials,100*mean(data.C==data.C_hat));