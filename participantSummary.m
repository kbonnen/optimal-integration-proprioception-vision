close all;
clear all;
choice = 'MLE';
nrStim = 2;

files = dir('./MLEcomb*');
id = 3434;
figure(id); clf;

for ii=1:15
    [results, sig,weights]=performIndividualAnalysisKLB(files(ii).name,[files(ii).folder '/'],choice,nrStim);
    
    
    
    % figure(1);
    % set(gcf, 'PaperSize', [4,8]);
    % set(gcf, 'PaperPosition', [0 0 4 8]);
    % saveas(gcf,[ 'individual/psychometric_' files(ii).name(1:end-4) '.pdf']);
    %
    % figure(2);
    % set(gcf, 'PaperSize', [8,4]);
    % set(gcf, 'PaperPosition', [0 0 8,4]);
    % saveas(gcf,[ 'individual/sigma_' files(ii).name(1:end-4) '.pdf']);
    
    figure(id);
    color = [];
    switch  char(sig + 48)
        case '111'
            text = 'optimal';
            color = 'k.';
            
        case '101'
            text = 'near-optimal';
            color = 'b.';
            
        case '011'
            text = 'ambiguous';
            color = 'c.';
            
        case '110'
            color = 'k.';
            
        case '001'
            text = 'not optimal';
            color = 'r.';
            
    end
    
    subplot(131); hold on;
    plot(results.cueCombos{1}.predicted_weights(1),results.cueCombos{1}.weights(1), color,'Markersize',20)
    xlim([0,1.5]); ylim([0,1.5]);
    axis square;
    
    subplot(132); hold on;
    plot(results.cueCombos{1}.predicted_sigma,results.cueCombos{1}.sigma, color,'Markersize',20)
    hold on;
    xlim([.5,4]); ylim([.5,4]);
    axis square;
    
    subplot(133); hold on;
    plot(results.singleCues{1}.params(2),results.singleCues{2}.params(2),color,'MarkerSize',20)
    hold on;
    xlim([0,5.5]); ylim([0,5.5]);
    xlabel('vision');
    ylabel('proprioception');
    axis square;
    
    close(figure(1));
    close(figure(2));
    
    
end

%%
subplot(131); refline(1,0);
title('weights');
xlabel('predicted');
ylabel('actual');

subplot(132); refline(1,0);
title('sigma');
xlabel('predicted');
ylabel('actual');

subplot(133); refline(1,0);

