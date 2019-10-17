% function entropy_simulation_fig_5A
% 
% Generates Figure 5A from "Optimal inference of sameness" by Van den Berg, 
% Vogel, Josic, and Ma, PNAS 2012.

% Written by RvdB, Feb 2012

function entropy_simulation_fig_5A

% parameters
sigma_s   = 5;
sigma_int = 10;
N=16;

% define partitions (based on Table 2 from Young & Wasserman 1997, Journal of Experimental Psychology)
parts{1} = [8,8];
parts{2} = [4,4,4,4];
parts{3} = [2,2,2,2,2,2,2,2];
parts{4} = [1 1 14];
parts{5} = [1 1 1 1 12];
parts{6} = [1 1 1 1 1 1 1 1 8];
parts{7} = [1 1 1 1 1 1 1 1 1 1 1 1 4];
parts{8} = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 2];
parts{9} = [13 1 1 1];
parts{10} = [10 3 2 1];
parts{11} = [7 5 3 1];
parts{12} = [4 4 4 4];
parts{13} = ones(1,16);
parts{14} = [16];

nTrials = 1000;  % number of trials per entropy value

% compute decision criterion (this is independent of the partitioning)
d = (sigma_int^2/N) * (1+sigma_int^2/sigma_s^2) * ((N-1) * log(1+sigma_s^2/sigma_int^2));
E = []; R = [];

for ii=1:length(parts)
    Q = parts{ii};
    for jj=1:nTrials
        % create stimulus set
        S = [];
        for kk=1:length(Q)           
            randOrt = round(randn*sigma_s); % draw an orientation from N(0,sigma_s) and round to nearest integer
            while ~isempty(find(S==randOrt, 1)) % redraw if any of the other subsets already has this orientation, because subsets should be distinct
                randOrt = round(randn*sigma_s);
            end
            S = [S ones(1,Q(kk))*randOrt]; % add subset to stimulus set
        end
        
        % compute entropy (in same way as Young & Wasserman did) and the ideal observer's response on this trial 
        E(end+1) = -sum(Q/N .* log2(Q/N));  % entropy of stimulus set
        x = S+randn(1,N)*sigma_int;         % simulate internal representation
        R(end+1) = var(x,1)>d;              % compute response
        
    end
end

uE = unique(E);     % determine set of unique entropy values 
uE = uE(2:end-1);   % we don't compute y-value for first and last bin ("all same" and "all different")

% compute "scaled logit of percent reported different" (as done by Young and Wasserman)
p_same = mean(R(E==0));  % proportion of "different" responses on trials where all items are the same
p_diff = mean(R(E==max(E(:)))); % proportion of "different" responses on trials where all items are different
for ii=1:length(uE)
    p_target=mean(R(E==uE(ii)));
    p(ii) = (p_target-p_same)/(p_diff-p_same);
    cnt(ii) = sum(E==uE(ii));
end

%  plot results and linear fit
close all
figure
set(gcf,'Position',get(gcf,'Position').*[.1 .1 .7 1]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 .7 1]);
set(gca,'FontSize',12,'FontName','Arial');
hold on;
plot(uE,p,'kx','MarkerSize',9,'LineWidth',2);
plot([0 10],polyval(polyfit(uE,p,1),[0 10]),'k');  % plot linear fit
xlim([min(uE)*.9 max(uE)*1.1])
xlabel('Entropy of Mixture Array');
ylabel('Scaled logit of percent different responses');
axis([0 4 0 1.1]);
