%% Comments

% Run this code to generate results that are then used to make tables and
% figures in the paper.

clc
clear
close
prompt = "This code will run all the experiments and then store the results. In order to see the results in the paper, run " + ...
    "the 'showResults.m'. Otherwise, press 'y' to run the whole simulations and store the results. \n";
answer = input(prompt, "s");

%% Main code

if strcmp(answer, 'y') || strcmp(answer, 'Y') || strcmp(answer, 'yes') || strcmp(answer, 'Yes')
    pause(3);
    clc
    close
    clear
    baseDir = pwd;
    noise_types = {'betaPrime', 'bimodal', 'cauchy', 'exponential', 'gamma', 'impulsive', 'levy', 'skewedNormal'};
    numMeasNoise = length(noise_types);
    for iii=1:numMeasNoise
        pause(3);
        clc
        close
        clearvars -except iii numMeasNoise baseDir noise_types
        inner = noise_types{iii};
        mainFile = fullfile(baseDir, inner, 'main.m');
        if exist(mainFile, 'file')
            fprintf('Running: %s\n', mainFile);
            try
                run(mainFile);
            catch ME
                warning('Failed: %s\n  -> %s', mainFile, ME.message);
            end
        else
            warning('Missing file: %s', mainFile);
        end
    end
end
clc