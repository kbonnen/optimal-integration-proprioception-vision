clear all;

inits = 'A';  % Change this for the participants (A, MW)
visual_noise_level = -1;


filenames = {['MLE_' inits '35.mat'], ['MLE_' inits '65.mat'], ['MLE_' inits '_S2.mat']};
trials = [];
for ff=1:3
    load(filenames{ff});
    trials = [trials; MLE];
end
trials(:,11) = [35*ones(156,1);65*ones(156,1);-1*ones(312,1)];


% find the unimodal trials
vision_ind = trials(:,5)==1 & trials(:,6)==0 & trials(:,11)==visual_noise_level; 
vision_tmp = (trials(vision_ind,:));
vision = table(vision_tmp(:,1),vision_tmp(:,2),vision_tmp(:,11),'VariableNames',{'level','choice','noise'});
unimodal = {};
g = grpstats(vision,{'level'},{'sum','mean'});
unimodal{1} = [g.level,g.sum_choice,g.GroupCount];

proprioceptive_ind = trials(:,6)==1 & trials(:,5)==0;
proprioceptive_tmp = (trials(proprioceptive_ind,:));
proprioceptive = table(proprioceptive_tmp(:,1),proprioceptive_tmp(:,2),'VariableNames',{'level','choice'});
g = grpstats(proprioceptive,{'level'},{'sum','mean'});
unimodal{2} = [g.level,g.sum_choice,g.GroupCount];



%% use psignifit to estimate sigma values

conditions = {'vision', 'proprioception'};
figure(2); clf;
for ii=1:2
    subplot(2,2,ii); hold on;
    psigni_data = unimodal{ii};
    options             = struct;   % initialize as an empty struct
    options.sigmoidName = 'norm';
    options.expType     = 'YesNo';
    options.estimateType = 'MAP';  % THis actually seems to matter more than I'd like.  MAP vs mean  estimator
    options.fixedPars = [nan,nan,nan,nan,0]';
    
    % Now run psignifit
    result{ii} = psignifit(psigni_data,options);
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
    plotOptions.fontSize       = 16;
    plotOptions.labelSize      = 16;
    
    plotPsych(result{ii},plotOptions);
    
    p = getStandardParameters(result{ii});
    
    sig(ii) = p(2);
    PSE(ii) = p(1);
    title(conditions{ii})
  
%     jnd(ii) = getThreshold(result,.84)-getThreshold(result,.5); %sqrt(2*sig(ii).^2);
end

% adjust the visual sigma slightly up/down ?
% sig(1) = 2


sig_bi_predicted = sqrt(prod(sig.^2)/sum(sig.^2));
sig(3) = sig_bi_predicted;  
PSE(3) = 0; % update this with the predicted code...

subplot(224);
for ii=1:2
   fprintf('%s: PSE=%.02f, sig=%.02f\n',conditions{ii},PSE(ii),sig(ii)); 
end
fprintf('predicted bimodal sigma value: %.02f\n',sig_bi_predicted);

plot(sig(1:2),'k.-','MarkerSize',20)
ylim([1,4]);
xlim([.5,2.5]);
set(gca','XTick',[1,2],'XTickLabels',{'v','p'});
xlabel('visual noise');
ylabel('sigma');


%% Ask the question of whether the smallest unimodal sigma and 
%  the bimodal prediction are distinguishable



levels = g.level;
dist = ones(size(levels))/numel(levels);
nTrials = [468+78,486+78,468*2];
nRepeats = 100;
figure(3); clf;

for ii=1:3
    ii
% now simulate from the mean predicted sigma to see what the error bar on
% the data would be.
    trialCount = round(nTrials(ii)*dist);
    for rr=1:nRepeats
        options.fixedPars = [nan,nan,0,0,0]';
        p = normcdf(levels,PSE(ii),sig(ii));
        choice = binornd(trialCount,p);
        data = [levels,choice,trialCount];
        result = psignifit(data,options);
        params(rr,:,ii) = result.Fit;
        params_CI(rr,:,:,ii) = result.conf_Intervals(:,:,1);
        s = getStandardParameters(result);
        bootstrap_sig(rr,ii) = s(2);
    end
    subplot(121); hold on;
    
    plotOptions.dataColor      = colors(ii,:);  % color of the datapoints
    plotOptions.plotData       = true;
    plotOptions.lineColor      = colors(ii,:);          % Color of the psychometric function
    plotOptions.lineWidth      = 2;                % Thikness of the psychometric function
    plotOptions.CIthresh       = true;            % plot a confidence interval at threshold
    plotOptions.fontSize       = 16;
    plotOptions.labelSize      = 16;
    plotPsych(result,plotOptions);

end

text(-8,.5,num2str(sig','%.02f'));


subplot(122); hold on;
xx=1:3;
yy = mean(bootstrap_sig);
neg = quantile(bootstrap_sig,.025);
pos = quantile(bootstrap_sig,.975);

for ii=1:3
plot(xx(ii),yy(ii),'.','MarkerSize',20,'Color',colors(ii,:));
e = errorbar(xx(ii),yy(ii),yy(ii)-neg(ii),pos(ii)-yy(ii));
set(e,'MarkerSize',20,'LineWidth',2,'Color',colors(ii,:));
end
xlim([.5,3.5]);

