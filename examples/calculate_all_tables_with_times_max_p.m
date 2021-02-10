%% All table generation script with computation times (max. prevalence)
% 
% This script uses the Toxo library to calculate penetrance tables of all
% the models available in this repository (i.e. additive, multiplicative 
% and threshold, from 2 to 8 order), using MAFs 0.1 and 0.4 for all locus 
% and heritabilities 0.1 and 0.8. This version maximizes prevalence.
%
% Warning: some models could not be computed to the given MAFs and
% heritabilities combination. Warnings are raised to explicit this.
%
% Author: Borja González Seoane
% Last modified: 9 Feb 2021
% Based in Christian Ponte's 'calculate_tables.m' script.

%% Script
% Import Toxo library
addpath('src/');

% Adjust Matlab accuracy. Necessary to compute some models without errors
digits(64);

% Read all the models from 'models' folder
model_list = {
    'models/additive_3.csv'
    'models/additive_3.csv'
    'models/additive_4.csv'
    'models/additive_5.csv'
    'models/additive_6.csv'
    'models/additive_7.csv'
    'models/additive_8.csv'
    'models/multiplicative_2.csv'
    'models/multiplicative_3.csv'
    'models/multiplicative_4.csv'
    'models/multiplicative_5.csv'
    'models/multiplicative_6.csv'
    'models/multiplicative_7.csv'
    'models/multiplicative_8.csv'
    'models/threshold_2.csv'
    'models/threshold_3.csv'
    'models/threshold_4.csv'
    'models/threshold_5.csv'
    'models/threshold_6.csv'
    'models/threshold_7.csv'
    'models/threshold_8.csv'
    };

% Create a list of MAFs and heritabilities to test
maf = [0.1, 0.4];
h = [0.1, 0.8];

% Create the output folder
output_folder = "output";
if ~ isfolder(output_folder)
    mkdir(output_folder);
end

% Main
for i = 1:length(model_list)
    % Create a Model instance
    m = toxo.Model(model_list{i});
    for j = 1:length(maf)
        mafs = repmat(maf(j), 1, m.order);
        for k = 1:length(h)
            try
                tic;  % To time the table calculation
                % Calculate a penetrance table for the model
                pt = m.find_max_prevalence(mafs, h(k));
                comptime = toc;  % To time the table calculation
                % Write the table to a file
                file_name = sprintf("%s_%.1f_h%.1f.csv", m.name, maf(j), h(k));
                pt.write(fullfile(output_folder, file_name), toxo.PTable.format_csv);
                % Save the computation time into the times file
                fid = fopen(fullfile(output_folder, 'times.txt'), 'a+');
                fprintf(fid, "%s_%.1f_h%.1f: %f\n", m.name, maf(j), h(k), comptime);
                fclose(fid);
            catch ME
                warning("(%s with MAF=%.1f and h²=%.1f) %s", m.name, maf(j), h(k), ME.message);
                continue;
            end
        end
    end
end
