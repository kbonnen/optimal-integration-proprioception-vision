close all;
clear all;
pathname = '/Users/kathrynbonnen/Documents/work-repos/TrevorData2022/';
filename = 'MLEcomb1.mat';
choice = 'MLE';
nrStim = 2;

%%
[results, ~,~]=performIndividualAnalysisKLB(filename,pathname,'MLE',nrStim);

%%
load(filename);
unimodal = {}
% find the unimodal trials
vision_ind = MLEcomb(:,5)==1 & MLEcomb(:,6)==0;
vision_tmp = (MLEcomb(vision_ind,:));
unimodal{1} = table(vision_tmp(:,1),vision_tmp(:,2),'VariableNames',{'level','choice'});


proprioceptive_ind = MLEcomb(:,6)==1 & MLEcomb(:,5)==0;
proprioceptive_tmp = (MLEcomb(proprioceptive_ind,:));
unimodal{2} = table(proprioceptive_tmp(:,1),proprioceptive_tmp(:,2),'VariableNames',{'level','choice'});




%% Bootstrap analysis of unimodal data to get sigma values and confidence intervals
% figure(1); subplot(211); hold on;
nRepeats = 220;
sig = nan(nRepeats,2);
PSE = nan(nRepeats,2);

for jj=1:nRepeats
    fprintf(num2str(jj))
    tic;
    ind = randsample(220,220,1);
    for ii=1:2
        
        % choose 100 with replacement
        
        raw = unimodal{ii}(ind,:);
        g = grpstats(raw,{'level'},{'sum'});
        
        psigni_data = [g.level, g.sum_choice, g.GroupCount];
        options             = struct;   % initialize as an empty struct
        options.sigmoidName = 'norm';
        options.expType     = 'YesNo';
        options.estimateType = 'mean';  % THis actually seems to matter more than I'd like.  MAP vs mean  estimator
%         lapse = [g.sum_choice(1)/g.GroupCount(1), 1-g.sum_choice(end)/g.GroupCount(end)];
        options.fixedPars = [nan,nan,0,0,0]';
        
        % Now run psignifit
        result = psignifit(psigni_data,options);
        p = getStandardParameters(result);
        
        sig(jj,ii) = p(2);
        PSE(jj,ii) = p(1);
        
    end
    toc;
end

%% Summary of estimated unimodal sigmas and predicted multi-modal sigma

figure(2); clf;
subplot(121);
histogram(sig(:,1),0:.125:6); hold on;
histogram(sig(:,2),0:.125:6);
sig_predicted = sqrt((sig(:,1).^2.*sig(:,2).^2)./sum(sig.^2,2));
histogram(sig_predicted,0:.125:6);
title('sigma')
legend({'vision','proprioception','predicted - combo'});

subplot(122); hold on;
mean(sig)
mean(sig_predicted)
% plot the actual
xx = [1,2];
yy = mean(sig);
neg = quantile(sig,.025);
pos = quantile(sig,.975);
plot(xx,yy,'k.','MarkerSize',20);
errorbar(xx,yy,yy-neg,pos-yy,'k.')
set(gca,'XTick',[1,2,3],'XTickLabels',{'v','p','b'});
xlim([0,4]);
title('sigma - summary');

yy = mean(sig_predicted);
neg = quantile(sig_predicted,.025);
pos = quantile(sig_predicted,.975);
plot(3,yy,'r.','MarkerSize',20);
errorbar(3,yy,yy-neg,pos-yy,'r')

%% Let's do a power analysis for this single subject and the predicted cue combination

% Number of trials to simulate
N = 440;


% set the levels for the experiment
levels = [-10,-5:5,10]';
% levels = [-10, -8, -6, -5, -4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10]';
% levels = [-10:2:10]';

dist = ones(size(levels)); %uniform
% dist = [10,1,2,4,8,16*ones(1,9),8,4,2,1,10]'; %increasing

% g = grpstats(unimodal{1},{'level'},{'sum'});
% dist = g.GroupCount; % approximately the staircase that someone experienced
% levels = g.level;

dist = dist/sum(dist);
trialCount = round(N*dist);

nRepeats = 200; % number of times to repeat the simulation
simulated_sig_unimodal = nan(nRepeats,2);

for rr = 1:nRepeats
    rr
    options.fixedPars = [nan,nan,nan,nan,0]';
    for jj=1:2
        p = normcdf(levels,mean(PSE(:,jj)),mean(sig(:,jj)));
        choice = binornd(trialCount,p);
        data = [levels,choice,trialCount];
        result = psignifit(data,options);
        s = getStandardParameters(result);
        simulated_sig_unimodal(rr,jj) = s(2);
    end
end

q = quantile(simulated_sig_unimodal,[.05,.5,.95]);
simulated_sig_predicted = sqrt((simulated_sig_unimodal(:,1).^2.*simulated_sig_unimodal(:,2).^2)./sum(simulated_sig_unimodal.^2,2));

% now simulate from the mean predicted sigma to see what the error bar on
% the data would be.
combo_sig = nan(nRepeats,1);
trialCount = round(2*N*dist);
for rr=1:nRepeats
    options.fixedPars = [nan,nan,0,0,0]';
    p = normcdf(levels,0,mean(simulated_sig_predicted));
    choice = binornd(trialCount,p);
    data = [levels,choice,trialCount];
    result = psignifit(data,options);
    params(rr,:) = result.Fit;
    params_CI(rr,:,:) = result.conf_Intervals(:,:,1);
    s = getStandardParameters(result);
    combo_sig(rr) = s(2);
    close all;
end



figure(3); clf;
hold on;
% plot the actual
xx = [1,2];
yy = mean(simulated_sig_unimodal);
neg = quantile(simulated_sig_unimodal,.025);
pos = quantile(simulated_sig_unimodal,.975);
plot(xx,yy,'k.','MarkerSize',20);
errorbar(xx,yy,yy-neg,pos-yy,'k.')
set(gca,'XTick',[1,2,3],'XTickLabels',{'v','p','b'});
xlim([0,4]);

yy = mean(simulated_sig_predicted);
neg = quantile(simulated_sig_predicted,.025);
pos = quantile(simulated_sig_predicted,.975);
plot(3,yy,'r.','MarkerSize',20);
errorbar(3,yy,yy-neg,pos-yy,'r')

yy = mean(combo_sig);
neg = quantile(combo_sig,.025);
pos = quantile(combo_sig,.975);
plot(3,yy,'g.','MarkerSize',20);
errorbar(3,yy,yy-neg,pos-yy,'g')




title('sigma - summary');
