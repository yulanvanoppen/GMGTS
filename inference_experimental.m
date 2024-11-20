%% Load data ---------------------------------------------------------------
clearvars
close all
rng('default');
warning('off','MATLAB:table:ModifiedAndSavedVarnames')


data = struct;


data(1).file = 'FP_data/EL222_AQTrip-mTurquoise2_microfluidic_CFP_15x300ms_rep1.xlsx';
DATA = readtable(data(1).file);
data(1).y = reshape(DATA.CFP_total ./ DATA.volume, 91, []);
data(1).kdil = .0033;

data(2).file = 'FP_data/EL222_AQTrip_mTFP1_microfluidic_GFP(new)_12x300ms_rep1.xlsx';
DATA = readtable(data(2).file);
DATA = [DATA(1:500, :); DATA(583:end, :)];
data(2).y = reshape(DATA.GFP_total ./ DATA.volume, 100, []);
data(2).kdil = .0038;

data(3).file = 'FP_data/EL222_AQTrip-mScarletI_microfluidic_CFP_1x1s_rep1.xlsx';
DATA = readtable(data(3).file);
data(3).y = reshape(DATA.RFP_total ./ DATA.volume, 120, []);
data(3).kdil = .0043;

data(4).file = 'FP_data/EL222_AQTrip-mCherry_microfluidic_CFP_1x1s_rep2.xlsx';
DATA = readtable(data(4).file);
DATA = [DATA(1:6, :); DATA(8:10, :); DATA(12:end, :)];
data(4).y = reshape(DATA.RFP_total ./ DATA.volume, 120, []);
data(4).kdil = .0048;

data(5).file = 'FP_data/EL222_AQTrip_tdTomato_microfluidic_CFP_1x1s_rep1';
DATA = readtable(data(5).file);
data(5).y = reshape(DATA.RFP_total ./ DATA.volume, 81, []);
data(5).kdil = .0042;

data(6).file = 'FP_data/EL222_AQTrip-mKate2_microfluidic_CFPnew_1x1s_rep2.xlsx';
DATA = readtable(data(6).file);
DATA = DATA(DATA.TimeID <= 120, :);
DATA = [DATA(1:3600, :); DATA(3720:end, :)];
data(6).y = reshape(DATA.RFP_total ./ DATA.volume, 120, []);
data(6).kdil = .0044;



data(7).file = 'FP_data/EL222_AQTrip-CFP_microfluidic_CFP_15x300ms_rep1.xlsx';
DATA = readtable(data(7).file);
data(7).y = reshape(DATA.CFP_total ./ DATA.volume, 90, []);
data(7).kdil = .0038;

data(8).file = 'FP_data/EL222_AQTrip-sfGFP_microfluidic_CFP_1x1s_rep1.xlsx';
DATA = readtable(data(8).file);
data(8).y = reshape(DATA.GFP_total ./ DATA.volume, 120, []);
data(8).y = data(8).y(:, ~ismember(1:size(data(8).y, 2), [3 12 13 19 20 22 42]));
data(8).kdil = .0042;

data(9).file = 'FP_data/EL222_AQTrip-pHtdGFP_microfluidic_GFP(new)_12x300ms_rep1.xlsx';
DATA = readtable(data(9).file);
data(9).y = reshape(DATA.GFP_total ./ DATA.volume, 90, []);
data(9).kdil = .0038;

data(10).file = 'FP_data/EL222_AQTrip-mVenus_microfluidic_CFP_1x1s_rep2.xlsx';
DATA = readtable(data(10).file);
DATA = DATA([1:712 795:end], :);
data(10).y = reshape(DATA.YFP_total ./ DATA.volume, 89, []);
data(10).kdil = .0045;

data(11).file = 'FP_data/EL222_AQTrip-mCitrine_microfluidic_CFP_1x1s.xlsx';
DATA = readtable(data(11).file);
DATA = DATA(DATA.TimeID <= 90, :);
DATA = [DATA(1:630, :); DATA(897:end, :)]; 
data(11).y = reshape(DATA.YFP_total ./ DATA.volume, 90, []);
data(11).kdil = .0043;

data(12).file = 'FP_data/EL222_AQTrip-mNeonGreen_microfluidic_CFP_1x1s_rep1.xlsx';
DATA = readtable(data(12).file);
data(12).y = reshape(DATA.YFP_total ./ DATA.volume, 120, []);
data(12).kdil = .0041;



for idx = 1:length(data)
    data(idx).name = regexprep(data(idx).file, "FP_data/EL222_AQTrip[-_]", "");
    data(idx).name = regexprep(data(idx).name, "_.+", "");
    t = 0:5:210;
    ind = 1:23+20*(idx<7);
    y = interp1(t(ind)', data(idx).y(ind, :), t(ind(1:end-2))+6);
    y = max(0, y - y(1, :));
    data(idx).y = y;
    data(idx).t = t(ind(1:end-2));
    data(idx).T = length(data(idx).t);
    data(idx).N = size(y, 2);
    data(idx).traces = reshape(y, size(y, 1), 1, size(y, 2));
end



%% Estimate ----------------------------------------------------------------
methods = [];
methods = [methods "GMGTS"];
% methods = [methods "GTS"];

knots = [12.5 25 50];

mkdir('estimates_experimental')
mkdir('estimates_experimental_fixed_km')

for idx = size(DATA, 2)
    if idx < 7, model = 'model_maturation_twostep.txt'; else, model = 'model_maturation_onestep.txt'; end
    disp(data(idx).name);
    
    system = System(model, 'FixedParameters', ["kr" "kdr" "kdil" "d"]);
    system.fixed.values(3) = data(idx).kdil;
    data(idx).init = system.x0' + 1e-6;
    data(idx).observed = system.K;
    estimator = Estimator(system, data(idx) ...                             % estimator setup
                          , 'Stages', 2 ...                                 % 0: smoothing only, 1: first stage only
                          , 'Methods', methods ...                          % GMGT, GTS, or both
                          , 'Knots', knots ...
                          , 'LB', [.001 .001] ...
                          , 'UB', [20 1] ...
                          , 'LogNormal', true ...
                         );
    estimator.estimate();
%     data(idx).logL_variable_GMGTS = estimator.loglik(5);
%     data(idx).logL_variable = estimator.loglik(5, "GTS");
    
    close all
    plot(estimator, 'States', 1:system.K, 'MaxCells', 100)
    
% % % %     save(['estimates_experimental/'  data(idx).file(9:end-5) '.mat'])
    

% % % % % % %     system = System(model, 'FixedParameters', ["kr" "kdr" "kdil" "d" "km"]);
% % % % % % %     system.fixed.values(end) = estimator.results_GMGTS.b_est(end);
% % % % % % %     estimator = Estimator(system, data(idx) ...                             % estimator setup
% % % % % % %                           , 'Stages', 2 ...                                 % 0: smoothing only, 1: first stage only
% % % % % % %                           , 'Methods', methods ...                          % GMGT, GTS, or both
% % % % % % %                           , 'Knots', knots ...
% % % % % % %                           , 'LB', .001 ...
% % % % % % %                           , 'UB', 20 ...
% % % % % % %                          );
% % % % % % %     estimator.estimate();
% % % % % % % %     data(idx).logL_fixed_GMGTS = estimator.loglik(5);
% % % % % % % %     data(idx).logL_fixed = estimator.loglik(5, "GTS");
% % % % % % %     
% % % % % % % %     plot(estimator, 'States', 2:system.K, 'MaxCells', 10)
% % % % % % %     
% % % % % % % % % % %     save(['estimates_experimental_fixed_km/'  data(idx).file(9:end-5) '.mat'])
end



%% Likelihood comparison
% [std(reshape([data.logL_fixed_GMGTS], [], 12)', 0, 2),...
%  mean(reshape([data.logL_fixed_GMGTS], [], 12)', 2),...
%  mean(reshape([data.logL_variable_GMGTS], [], 12)', 2),...
%  std(reshape([data.logL_variable_GMGTS], [], 12)', 0, 2)]
% 
% mean(reshape([data.logL_fixed_GMGTS], [], 12)', 2) ...
%      < mean(reshape([data.logL_variable_GMGTS], [], 12)', 2)
%  
% printable = {data(:).name}';
% printable = [printable compose('%.0f', mean(reshape([data.logL_fixed_GMGTS], [], 12)', 2))];
% printable = [printable compose('%.0f', mean(reshape([data.logL_variable_GMGTS], [], 12)', 2))];
% fprintf(repmat('%s & %s & %s \\\\ \n', 1, length(data)), string(printable'))


