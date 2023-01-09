close all;
clear all;
pathname = '/Users/kathrynbonnen/Documents/work-repos/TrevorData2022/';
filename = 'MLEcomb2.mat';
choice = 'MLE';
nrStim = 2;

[results, sig,weights]=performIndividualAnalysisKLB(filename,pathname,choice,nrStim);

load('MLEcomb1.mat');

% find the unimodal trials
vision_ind = MLEcomb(:,5)==1 & MLEcomb(:,6)==0;
vision_tmp = (MLEcomb(vision_ind,:));
vision = table(vision_tmp(:,1),vision_tmp(:,2),'VariableNames',{'level','choice'});
g = grpstats(vision,{'level'},{'sum','mean'});
unimodal = {};
unimodal{1} = [g.level,g.sum_choice,g.GroupCount];

proprioceptive_ind = MLEcomb(:,6)==1 & MLEcomb(:,5)==0;
proprioceptive_tmp = (MLEcomb(proprioceptive_ind,:));
proprioceptive = table(proprioceptive_tmp(:,1),proprioceptive_tmp(:,2),'VariableNames',{'level','choice'});
g = grpstats(proprioceptive,{'level'},{'sum','mean'});
unimodal{2} = [g.level,g.sum_choice,g.GroupCount];



%%
% use psignifit to check  out the raw data

figure(1); subplot(211); hold on;
for ii=1:2
    
    psigni_data = unimodal{ii};
    options             = struct;   % initialize as an empty struct
    options.sigmoidName = 'norm';
    options.expType     = 'YesNo';
    options.estimateType = 'mean';  % THis actually seems to matter more than I'd like.  MAP vs mean  estimator
    options.fixedPars = [nan,nan,nan,nan,0]';
    
    % Now run psignifit
    result = psignifit(psigni_data,options);
    % [thresh0, thresh_ci0] = getThreshold(result,.5);
    % [slope0] = getSlope(result,thresh0);
    
    colors = [0, 0.4470, 0.7410; ...
        0.8500, 0.3250, 0.0980; ...
        0.9290, 0.6940, 0.1250];
    linetype = {'-','--'};
    
    plotOptions.dataColor      = colors(ii,:);  % color of the datapoints
    plotOptions.plotData       = true;
    plotOptions.lineColor      = colors(ii,:);          % Color of the psychometric function
    plotOptions.lineWidth      = 2;                % Thikness of the psychometric function
    plotOptions.CIthresh       = true;            % plot a confidence interval at threshold
    plotOptions.fontSize       = 32;
    plotOptions.labelSize      = 40;
    
    plotPsych(result,plotOptions);
    
    p = getStandardParameters(result);
    
    sig(ii) = p(2);
    PSE(ii) = p(1);
%     jnd(ii) = getThreshold(result,.84)-getThreshold(result,.5); %sqrt(2*sig(ii).^2);
end

sig_bi_predicted = sqrt(prod(sig.^2)/sum(sig.^2));