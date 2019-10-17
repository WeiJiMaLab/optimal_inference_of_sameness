% function blur_simulation_fig_5C
% 
% Generates Figure 5C from "Optimal inference of sameness" by Van den Berg, 
% Vogel, Josic, and Ma, PNAS 2012.

% Written by RvdB, Feb 2012

function blur_simulation_fig_5C

% init
sigma_s = 5;     % width of the distribution from which the "different" stimuli are drawn
sigma_int_vec = linspace(.5,70,19); % "blur" values (higher blur corresponds to lower sigma_int)
N = 8;           % set size
nTrials = 1000;  % number of trials to simulate per sigma_int value
psame = .6;      % prior for "same" displays
g = .25;         % guessing rate (i.e., proportion of trials on which subject makes random guess)

% simulate
for ii=1:length(sigma_int_vec)

    % compute decision criterion
    sigma_int = sigma_int_vec(ii);
    d = (sigma_int^2/N) * (1+sigma_int^2/sigma_s^2) * ((N-1)*log(1+sigma_s^2/sigma_int^2)+2*log(psame/(1-psame)));    
    
    % draw stimulus values 
    s_same = zeros(2,N);
    s_diff = normrnd(0,sigma_s,nTrials,N);

    % compute p(resp = "diff") for "same" and "diff" trials
    delta_vec_diff = sum(((s_diff-repmat(mean(s_diff,2),1,N)).^2 / sigma_int^2),2)';
    perf_diff(ii) = mean(1-ncx2cdf(N.*d/sigma_int^2,N-1,delta_vec_diff'));
    ncx2cdf(N.*d/sigma_int^2,N-1,delta_vec_diff');
    delta_vec_same = sum(((s_same-repmat(mean(s_same,2),1,N)).^2 / sigma_int^2),2)';
    perf_same(ii) = mean(ncx2cdf(N.*d/sigma_int^2,N-1,delta_vec_same'));       
    
end

% compute performance on same and on different trials
perf_same = (1-g)*perf_same + g/2;
perf_diff = (1-g)*perf_diff + g/2;

% plot result
figure;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 .75 1]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 .75 1]);
hold on;
plot(sigma_int_vec,100*perf_same,'ko-','MarkerSize',10);
plot(sigma_int_vec,100*perf_diff,'kx-','MarkerSize',10);
ylabel('Percent correct');
xlabel('\sigma');
set(gca,'YTick',0:25:100,'YTickLabel',{'0%','25%','50%','75%','100%'});
l=legend({'"Same" trials','"Different" trials'},'Location','East');
set(l,'box','off');
ylim([0 100])
box on;