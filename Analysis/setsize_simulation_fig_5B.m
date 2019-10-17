% function setsize_simulation_fig_5B
% 
% Generates Figure 5B from "Optimal inference of sameness" by Van den Berg, 
% Vogel, Josic, and Ma, PNAS 2012.

% Written by RvdB, Feb 2012

function setsize_simulation_fig_5B

% init
sigma_s = 8;                 % width of the distribution from which stimuli on "different" trials are drawn
sigma_int = 8;               % internal noise parameter
N_vec = [2 4 8 12 16 20 24]; % set sizes 
nTrials = 1000;              % number of trials to simulate per set size
psame = .47;                 % subject's prior for "same" displays

% simulate
for ii=1:length(N_vec)
    N = N_vec(ii);

    % compute decision criterion   
    d = (sigma_int^2/N) * (1+sigma_int^2/sigma_s^2) * ((N-1)*log(1+sigma_s^2/sigma_int^2)+2*log(psame/(1-psame)));    
    
    % draw stimulus values 
    s_same = zeros(2,N);
    s_diff = normrnd(0,sigma_s,nTrials,N);

    % compute p(resp = "diff") for "same" and "diff" trials
    delta_vec_diff = sum(((s_diff-repmat(mean(s_diff,2),1,N)).^2 / sigma_int^2),2)';
    p_diff_diff(ii) = mean(1-ncx2cdf(N.*d/sigma_int^2,N-1,delta_vec_diff'));
    ncx2cdf(N.*d/sigma_int^2,N-1,delta_vec_diff');
    delta_vec_same = sum(((s_same-repmat(mean(s_same,2),1,N)).^2 / sigma_int^2),2)';
    p_diff_same(ii) = mean(1-ncx2cdf(N.*d/sigma_int^2,N-1,delta_vec_same'));       
end

% plot
figure
N_vec(1) = 0; % just for plotting (plot N=2 at location of N=0, as in the Wasserman, Young, and Faggot paper)
plot(N_vec,100*p_diff_same,'ko-','MarkerSize',10);
hold on;
plot(N_vec,100*p_diff_diff,'kx-','MarkerSize',10);
hold off;
ylabel('Percent different responses');
xlabel('Number of stimuli');
set(gca,'YTick',0:25:100,'YTickLabel',{'0%','25%','50%','75%','100%'});
l=legend({'"Same" trials','"Different" trials'},'Location','East');
set(l,'box','off');
set(gca,'XTick',[0 4 8 12 16 20 24],'XTickLabel',[2 4 8 12 16 20 24]);
ylim([0 100])
xlim([-1 25])
box on;
 