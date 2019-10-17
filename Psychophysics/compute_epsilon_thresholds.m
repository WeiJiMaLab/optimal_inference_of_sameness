function [th_low th_high nBlocks] = compute_epsilon_thresholds(subjid,makeplot)

datadir = 'output';
files = dir(['output/' subjid '_condition_2*.mat']);

if isempty(files)
    th_low = -1;
    th_high = -1;
    nBlocks=-1;
    return;
end

allEpsilon = [];
allPerf = [];
nBlocks = length(files);
for ii=1:length(files)
    load([datadir '/' files(ii).name]);
    allEpsilon = [allEpsilon data(:,8)'];
    allPerf = [allPerf data(:,5)'];
end
nTrials = length(allPerf);

uEpsilon = unique(allEpsilon);

for ii=1:length(uEpsilon)
    X(ii) = uEpsilon(ii);
    idx   = allEpsilon == uEpsilon(ii);
    Y(ii) = sum(allPerf(idx)) / length(allPerf(idx));
end


% fit
fitpars = fminsearch(@(fitpars) fitfun(fitpars,uEpsilon,Y),[5 1 .4]);
fitpars(3) = min(.5,fitpars(3));    
X_fit = 0:.001:1;
Y_fit = 0.5 + fitpars(3) * normcdf(X_fit,fitpars(1),fitpars(2));
Y_fit = normrnd(Y_fit,.000001);   % to avoid "Values should be distinct" error in interp1

% low_prop =  0.5 + fitpars(3)*.50;
% high_prop = 0.5 + fitpars(3)*.99;

low_prop =  0.7;

th_low  = interp1(Y_fit,X_fit,low_prop); % set low to 70%-correct threshold
th_low  = nanmin(.84,th_low);   % set maximum low to .84
th_high = .94;               

if (makeplot)
    figure;
    set(gca,'FontSize',14);
    plot(X,Y,'rx','LineWidth',2,'MarkerSize',12);
    set(gca,'FontSize',12);
    hold on;
    plot(X_fit,Y_fit,'k-','LineWidth',2)
    legend({'data','cdf fit'},'Location','NorthWest');
    xlabel('\epsilon');
    text(.05,.9,['th_{low} = ' num2str(100*low_prop,2) '% threshold = ' num2str(th_low,2)]);
    text(.05,.87,['th_{high} = ' num2str(th_high,2) ' (fixed)']);
    text(.05,.84,['ceiling = ' num2str(.5+fitpars(3),2)]);
    ylabel('Proportion correct');
    ylim([min(min(Y),.5) 1.05]);
    xlim([0 1]);
    title(['Subject ' subjid ' (' num2str(nTrials) ' trials)']);
    plotfname = ['plots/epsilon_thresholds_' subjid '.png'];
    print(gcf,'-dpng',plotfname);   
end


% helper function
function sse = fitfun(pars,uEpsilon,Y)
pars(3)=min(.5,pars(3));
sse = sum( ((0.5 + pars(3) * normcdf(uEpsilon,pars(1),pars(2))) - Y).^2) ;