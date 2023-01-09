clear; close all;

levels{1} = [-10:10]';
dist{1} = ones(size(levels{1}));
dist{1} = dist{1}/sum(dist{1});

levels{2} = [-8,-6,-4:4,6,8]';
dist{2} =  ones(size(levels{2}));
dist{2} = dist{2}/sum(dist{2});

nRepeats = 100;
nTrials = 156;

figure(3); clf; 

for ii=2
% now simulate from the mean predicted sigma to see what the error bar on
% the data would be.
    sig = nan(nRepeats,1);
    trialCount = round(nTrials*dist{ii});
    for rr=1:nRepeats
        options.fixedPars = [nan,nan,0,0,0]';
        p = normcdf(levels{ii},0,2.5);
        choice = binornd(trialCount,p);
        data = [levels{ii},choice,trialCount];
        result = psignifit(data,options);
        params(rr,:) = result.Fit;
        params_CI(rr,:,:) = result.conf_Intervals(:,:,1);
        s = getStandardParameters(result);
        sig(rr) = s(2);
    end
    subplot(121);
    plotPsych(result)
    
    subplot(122);
    xx=ii
    yy = mean(sig)
    neg = quantile(sig,.025);
    pos = quantile(sig,.975);
    plot(xx,yy,'k.','MarkerSize',20);
    e = errorbar(xx,yy,yy-neg,pos-yy,'k.');
    set(e,'MarkerSize',20,'LineWidth',2);
    xlim([0,3]);

end