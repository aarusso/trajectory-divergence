%% DIVERGEANALYSIS.m
% This analysis determines whether a datasets possesses moments of high
% 'divergence'. Divergence will be considered high if there exist moments
% where the state of the system is similar but at some point later, the
% states become different. A constant (proportional to the variance of the
% data) is added to the denominator to prevent hyperbolic growth.
%
% Divergence is defined between every pair of times (t1 and t2) and is calculated
% for all times into the future (dT) until the end of the corresponding
% conditions. Because we are interested in whether datasets ever have these
% moments of high divergence, the value that is most informative is the
% maximum a given time point (t1) becomes divergent with any other time
% point (t2). We thus summarize divergence with a single value that
% indicates the maximum divergence for a given time point (t1) across all
% paired time points (t2) and dT. This value is the output 'D'.
%
% SYNTAX
%   [ D, params] = divergeAnalysis( Data ); 
%
% INPUTS:
%
% Required Inputs:
% * Data (struct) - a C-dimensional structure where C is the number of
% conditions. For a given condition, Data(c).A should hold the data (e.g.
% firing rates). Each column of A corresponds to a neuron, and each row to
% a timepoint. Data(c).times and Data(c).analyzeTimes are optional
% fields. size(Data(c).times,1) should equal size(Data(c).A,1). If
% Data(c).analyzeTimes is provided, only those entries that match between
% 'times' and 'analyzeTimes' will be used for the analysis. If analyzeTimes
% == [], all times will be used. 
% NB: data is formatted as for jPCA (Churchland et al 2012)
%
%
% Optional Inputs:
% Several parameters can be specified from the command line (see Usage
% Examples below)
% * numPCs (scalar, default: 12) - specifies how many PCs to use when
% calculating tangling
% * withinConditionsOnly (logical, default: false) - if desired, divergence
% will only be calculated for pairs of times that correspond to the same
% condition
% * timeStep (scalar, default: 20) - Because divergence is a pairwise
% computation, things can get quite slow if enough data is included. After
% calculating the derivitive, it can be efficient to calculate divergence
% only on every ith time sample where i is specified by timeStep. It is
% generally advisible to keep this number as low as possible without
% computation times getting unreasonably long. For the data in Russo et al
% 2020, computation time was < ~30 seconds and i was set to 20 (corresponding to
% 20ms). Note that if i is too large, aliasing could occur and divergence
% could appear lower than it actually is.
% * softenNorm (scalar, default: empty) - In the usual Churchlandian
% fashion, soft normalization is performed on the neural data such that
% each neuron is devided by its range (across all times and conditions) + a
% constant indicated by softenNorm. EMG data should be fully normalized
% (i.e. softenNorm set to 0). For neural data, softenNorm was set to 5
% which means that neurons with a firing rate of <5Hz will be somewhat
% discounted. In the Russo et al 2020 dataset, average firing rates were
% ~50Hz so softenNorm effectively fully normalized most neurons. Leave
% empty if data has already been normalized. If empty and data is not
% normalized, code will throw a warning
% * epsilon (scalar, default: .01) - This value corresponds to the constant
% added to the demonenator to prevent hyperbolic growth. That value will be
% epsilon * the variance of the data.
% * dTMax (scalar, default: Inf) - This value corresponds to how far into
% the future divergence will be considered. 
%
% OUTPUTS:
% * D (vector) - the divergence index for each time point,  a T/timesStep x 1 vector
% * params (struct) - out is a struct containing the settings applied to for
% calculation of divergence (epsion, numPCs, dMax, timeStep etc).
%
% USAGE EXAMPLES
%   [ D, params] = divergeAnalysis( Data, 'numPCs', 10 ); % to calculate divergence with 10 PCs
%   [ D, params] = divergeAnalysis( Data,'withinConditionsOnly', true,'numPCs', 12) % to calculate divergence with 12 PCs, only within conditions
%
% AUTHORSHIP
% |Author: Abigail Russo,|
% |Email: russoaa@gmail.com,|
% |Dated: May 2019|

function [ D, params ] = divergeAnalysis( Data, varargin )



params = inputParser;

% Set default options
params.addParameter('epsilon',  .01, @isscalar);
params.addParameter('dTMax', Inf, @isscalar);
params.addParameter('withinConditionsOnly', false, @islogical);
params.addParameter('timeStep', 20, @isscalar); % larger values will speed computation
params.addParameter('numPCs', 12, @isscalar);
params.addParameter('softenNorm', [], @(x) isvector(x) | isempty(x));

% Parse command line for parameter settings
params.parse(varargin{:});

numPCs = params.Results.numPCs;
epsilon = params.Results.epsilon;
withinConditionsOnly = params.Results.withinConditionsOnly;
timeStep = params.Results.timeStep;
dTMax = params.Results.dTMax;
softenNorm = params.Results.softenNorm;

D = [];

%% get data type of interest and mask for usable times

condIDNum = 1:length(Data);
A_cell = cell(size(condIDNum));
condMask = cell(size(condIDNum));


for cc = 1:length(Data)
   if ~isfield(Data(cc), 'times')
      Data(cc).times = 1:size(Data(cc).A,1);
   elseif isempty(Data(cc).times)
      Data(cc).times = 1:size(Data(cc).A,1);
   end
   
   if ~isfield(Data(cc), 'analyzeTimes')
      Data(cc).analyzeTimes = Data(cc).times;
   elseif isempty(Data(cc).analyzeTimes)
      Data(cc).analyzeTimes = Data(cc).times;
   end
end

% reformat and mask data for desired analysis times
for cc = condIDNum
   timeMask = ismember(Data(cc).times, Data(cc).analyzeTimes);
   A_cell{cc} = Data(cc).A(timeMask,:);
   condMask{cc} = cc*ones(size(Data(cc).analyzeTimes));
end


A = cell2mat(A_cell');

%% Normalize / mean center data
% soft-normalize firing rates in the usual Churchland way
if ~isempty(softenNorm)
   normFactors = range(A,1)+softenNorm;
   A = bsxfun(@times, A, 1./normFactors);
elseif any(range(A,1) > 1) && isempty(softenNorm)
   warning('A should be normalized or soft-normalized such that the range of each neuron <= 1')
end
A = bsxfun(@minus, A, mean(A,1)); % mean center


%% Perform dimensionality reduction
PCs = pca(A);

if numPCs > size(A,2)
   numPCs = size(A,2);
   disp('numPCs > number of recordings');
end
topPCs = PCs(:,1:numPCs);
X = A * topPCs;

% calculate normalization factor for the denomenator (proportional to
% variance of data)
alpha = epsilon*sum(var(X));
X = mat2cell(X, cellfun(@(s) size(s,1), A_cell), numPCs);

%% Compute divergence on each pair of conditions

numberOConds = length(condIDNum);

for cond1 = 1:numberOConds
   t1_l = [];
   t1_DT = [];
   
   for cond2 = 1:numberOConds
      if withinConditionsOnly && cond1~=cond2
         % if only computing divergence within conditions, don't compute
         % divergence for pairs of conditions that are not the same
      else
         [l, thisDT] = calcDivergence(X{cond1}, X{cond2});
         t1_l = [t1_l l]; 
         t1_DT = [t1_DT thisDT];
      end
   end
   
   % summarize divergence by taking the max across the second dimension
   [max_l, idx] = max(t1_l,[],2);
   D = [D max_l'];
end

%% Compute confusion
   function [divergence, best_dT] = calcDivergence(X1, X2)
      % calculate divergence for every pair of times between two conditions
      
      T1 = size(X1,1); T2 = size(X2,1);
      
      Future_state_diff = nan(T1,T2);
      Current_state_diff = nan(T1,T2);
      best_dT = nan(T1,T2);
      
      for t1 = 1:timeStep:T1
         for t2 = 1:timeStep:T2
            maxdT = min([T1 - t1, T2 - t2, dTMax]);
            
            Current_state_diff(t1,t2) = norm(X1(t1,:) - X2(t2,:),2)^2; % squared L2-norm
            
            future_diffs = X1(t1:t1+maxdT,:) - X2(t2:t2+maxdT,:);
            future_diffs = sqrt(sum((future_diffs).^2,2)); % the 2-norm
            future_diffs = future_diffs.^2; % squared L2-norm
            
            % take the max across dT
            [Future_state_diff(t1,t2), best_dT(t1,t2)] = max(future_diffs);
         end
      end
      
      Future_state_diff = Future_state_diff(1:timeStep:T1, 1:timeStep:T2);
      Current_state_diff = 1./(Current_state_diff(1:timeStep:T1, 1:timeStep:T2)+alpha);
      divergence = (Future_state_diff.*Current_state_diff);
      best_dT = best_dT(1:timeStep:T1, 1:timeStep:T2);
      
   end


%% save parameters

out = struct();

params2save = {'epsilon', 'withinConditionsOnly','timeStep', 'numPCs', 'dTMax'};
for pp = 1:length(params2save)
   out.(params2save{pp}) = eval(params2save{pp});
end

D = D';
end

