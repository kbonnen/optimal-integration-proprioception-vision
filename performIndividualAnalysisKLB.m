function[results, significance_sigma, significance_weight]=performIndividualAnalysisKLB(filename,pathname,choice,nrStim)
% clear all;
% close all;
% clc

% -------------------------------------------------------------------------
%
%   This file reads in data files from a cue integration experiment from a
%   directory and performs an analysis using the MLE_Analysis.m file.
%   please consult the documentation MLEtoolbox_documentation.pdf or type
%   'help MLE_Analysis' for information how to prepare the data.
%
%   Note that this analysis assumes you have used a 2AFC paradigm in which
%   a standard stimulus is compared to a comparison stimulus on each trial
%   (that is nrStim = 2).
%   Type "help MLE_Analysis" or "help MLE_SingleSubject" for further
%   options.
%
%   by Marieke Rohde, 2014
%   marieke.rohde@uni-bielefeld.de
%
%   patch added by Loes van Dam, 2015
%   - include nrStim question
%
% -------------------------------------------------------------------------


% query user for data directory

% [filename,pathname] = uigetfile('*.mat;*.txt','Please Select experimental results file');

% read in the file
try
    % search for files ending in .mat
    if sum(filename(end+(-3:0)) == '.mat') ==4
        data = load([pathname filename]);
        var = fieldnames(data);
        data = eval(['data.' var{1}]);
        
        % append data from all data files and add subject number as
        % additional column
        fprintf(['Reading in ' filename '\n\n']);
    end
    if sum(filename(end+(-3:0)) == '.txt') ==4
        data = load([pathname filename]);
        fprintf(['Reading in ' filename '\n\n']);
    end
catch
    disp(['There was a problem reading in the file ' filename ...
        ' - does it have the correct and identical format?'])
    return
end


% ask the user whether they want to just look at unimodal data
% or generate MLE predictions and bimodal analysis

% choice = questdlg('Just unimodal performance, or also MLE analysis?', ...
%     '', ...
%     'MLE','unimodal','MLE');


switch choice
    
    % user just wants to fit unimodal data
    case 'unimodal'
        fprintf('Fitting unimodal conditions for the specified data.\n\n')
        
        for cue = 1:2
            datacue = data((data(:,5)==mod(cue,2))& (data(:,6)==mod(cue+1,2)),:);
            noise = unique(datacue(:,8+cue));
            % generate one plot + fit for each noise level
            for i = 1:length(noise)
                fprintf(['cue ' int2str(cue) ' noise level ' num2str(noise(i)) '\n'])
                data_temp = datacue(datacue(:,8+cue)==noise(i),:);
                
                % fitting the data
                [mu,sigma,lapse,tmp] = FitCumulativeGauss(data_temp(:,1:3),0);
                
                % plotting the data
                comps = unique(data_temp(:,1));
                figure
                title(['cue ' int2str(cue) ' noise level ' num2str(noise(i))])
                hold on
                for j = 1:length(comps)
                    plot(comps(j),sum(data_temp(data_temp(:,1)==comps(j),2))./sum(data_temp(data_temp(:,1)==comps(j),3)),'.')
                end
                plot(linspace(comps(1),comps(end),200),normcdf(linspace(comps(1),comps(end),200),mu,sigma));
                xlabel('stimulus axis (standard-comparison)')
                ylabel('proportion comparison chosen')
                fprintf(['PSE: ' num2str(mu) ' JND: ' num2str(sigma) '\n\n'])
            end
        end
        % user wants to evaluate MLE cue combination data for a single subject
    case 'MLE'
     
        
    % ask the user how many stimuli were compared in each trial. E.g. if two
    % stimuli were compared the JND represents Sqrt(2)*noise in the signal.

%         nrStim = questdlg('Nr of stimuli/intervals compared in each trial? (usually 2)', ...
%     '', ...
%     '1','2','2');

        fprintf('Performing MLE analysis for the specified data.\n\n')
        
        % calculating results
        results = MLE_SingleSubject(data, (nrStim), 0,0,[1 1],1);
        
        % generate text output that summarizes the results on the screen
        
        fprintf('\n\n*************************** \n');
        fprintf('*** Results single cues *** \n');
        fprintf('*************************** \n\n');
        for i = 1:size(results.singleCues,2)
            for j = 1: size(results.singleCues{1,i}.noiselevels,2)
                fprintf(['Cue: ' int2str(i) ' noise level: ' num2str(results.singleCues{1,i}.noiselevels(j)) '\n\n'])
                fprintf([' PSE: ' num2str(results.singleCues{1,i}.params(j,1))  ',  conf. interval ['   ...
                    num2str(results.singleCues{1,i}.params_CI(j,1,1)) ' ' num2str(results.singleCues{1,i}.params_CI(j,1,2)) ']\n']);
                fprintf([' sigma: ' num2str(results.singleCues{1,i}.params(j,2))  ',  conf. interval ['   ...
                    num2str(results.singleCues{1,i}.params_CI(j,2,1)) ' ' num2str(results.singleCues{1,i}.params_CI(j,2,2)) ']\n\n']);
                sigmas(j,i,:)=results.singleCues{1,i}.params_CI(j,2,:) ;
            end
        end
        
        fprintf('\n****************************************** \n');
        fprintf('*** Predictions and results multimodal ***\n');
        fprintf('******************************************\n');
        for j = 1: size(results.cueCombos{1,1}.noiselevels,1)
            fprintf('\n** ')
            for i = 1: size(results.cueCombos{1,1}.noiselevels,2)
                fprintf(['noise cue ' int2str(i) ': ' num2str(results.cueCombos{1,1}.noiselevels(j,i)) '  '])
                
                sigmas_b(j,i,:) = sigmas((results.singleCues{1,i}.noiselevels == results.cueCombos{1,1}.noiselevels(j,i)),i,:);
            end
            fprintf('\n\n')
            
            
            fprintf('Weights: \n');
            for i = 1: size(results.cueCombos{1,1}.weights,2)
                fprintf(['\n cue ' int2str(i) ':\n   predicted ' num2str(results.cueCombos{1,1}.predicted_weights(j,i)) ...
                    ',  conf. interval [' num2str(results.cueCombos{1,1}.predicted_weights_CI(j,i,1)) ' ' num2str(results.cueCombos{1,1}.predicted_weights_CI(j,i,2)) ...
                    ']\n   empirical '  num2str(results.cueCombos{1,1}.weights(j,i)) ...
                    ',  conf. interval [' num2str(results.cueCombos{1,1}.weights_CI(j,i,1)) ' ' num2str(results.cueCombos{1,1}.weights_CI(j,i,2)) ']\n'])
                
                %% not significantly different from prediction?
                significance_weight(j,i)= results.cueCombos{1,1}.predicted_weights_CI(j,i,1)<results.cueCombos{1,1}.weights_CI(j,i,2) ...
                    && results.cueCombos{1,1}.predicted_weights_CI(j,i,2)>results.cueCombos{1,1}.weights_CI(j,i,1) ;
            end
            fprintf('\n');
            
            fprintf('Sigma: \n\n');
            fprintf(['   predicted ' num2str(results.cueCombos{1,1}.predicted_sigma(j))  ...
                ',  conf. interval [' num2str(results.cueCombos{1,1}.predicted_sigma_CI(j,1)) ' ' num2str(results.cueCombos{1,1}.predicted_sigma_CI(j,2)) ']']);
            fprintf(['\n   empirical ' num2str(results.cueCombos{1,1}.sigma(j))  ...
                ',  conf. interval [' num2str(results.cueCombos{1,1}.sigma_CI(j,1)) ' ' num2str(results.cueCombos{1,1}.sigma_CI(j,2)) ']\n'])
            fprintf('\n');
            
            %% significantly smaller than better unimodal?
            significance_sigma(j,1) = results.cueCombos{1,1}.sigma_CI(j,2)<min(sigmas_b(j,:,1),[],2);
            
            %% not significantly larger than prediction?
            significance_sigma(j,2) = results.cueCombos{1,1}.sigma_CI(j,1)<results.cueCombos{1,1}.predicted_sigma_CI(j,2);
            
            %% not significantly smaller than prediction?
            significance_sigma(j,3) = results.cueCombos{1,1}.sigma_CI(j,2)>results.cueCombos{1,1}.predicted_sigma_CI(j,1);
            
        end
        
        fprintf('\n***************************************** \n');
        fprintf('*** Significance of results (summary) ***\n');
        fprintf('*****************************************\n');
        for j = 1: size(results.cueCombos{1,1}.noiselevels,1)
            fprintf('\n** ')
            for i = 1: size(results.cueCombos{1,1}.noiselevels,2)
                fprintf(['noise cue ' int2str(i) ': ' num2str(results.cueCombos{1,1}.noiselevels(j,i)) '  '])
            end
            fprintf('\n')
            
            for i = 1: size(results.cueCombos{1,1}.weights,2)
                fprintf(['\n weigth cue ' int2str(i) ': '])
                if significance_weight(j,i)
                    fprintf(' optimal ');
                else
                    fprintf(' not optimal ');
                end
            end
            
            text = '';
            switch  char(significance_sigma(j,:) + 48)
                case '111'
                    text = 'optimal';
                case '101'
                    text = 'near-optimal';
                case '011'
                    text = 'ambiguous';
                case '110'
                    text = 'super-optimal';
                case '001'
                    text = 'not optimal';
            end
            
            fprintf(['\n\n sigma ' text])
            
            fprintf('\n\n')
            
        end
        fprintf('(see Rohde, van Dam & Ernst, Multisensory Research (2015) for explanation of the labels)\n\n');
       
end