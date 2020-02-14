function [W,H,measurementIndexTable] = buildMeasurementMatrices(obsTimes, ...
                                                    kalmanOptions, ...
                                                    numVar,measurements)
  
% function [W,H,measurementIndexTable] =
% buildMeasurementMatrices(obsTimes,kalmanOptions,...  numVar,measurements)
%
% Calculate measurement covariance matrix (W) and measurement index
% matrix (H), based on either given weights (in kalmanOptions) or based
% on tables. The tables consist of two vectors - measurementIndex and
% measurementLabel - and one set of tables - sigmaVectors. Entry i in
% sigmaVectors contains data from sensor given by the name in
% measurementLabel{i}. measurementIndex(i) corresponds to the number in
% which it is defined in the SUMMARY section in the main eclipse data
% file.   
%  
% Input:
% ------
% obsTimes      : observation times
% kalmanOptions : kalman options
% numVar        : total number of variables 
% measurements  : measurements at obsTimes
%
% Output:
% -------
% W                     : measurement covariance matrix
% H                     : measurement index matrix
% measurementIndexTable : table indicating which measurements each W and
%                         H consists of
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/buildMeasurementMatrices.m#1 $
% $DateTime: 2011/10/13 14:36:15 $
 

% number of measurements
numMeas = size(measurements,1);

% use table to build H and W
if isfield(kalmanOptions,'measurementTable') && ...
   kalmanOptions.measurementTable == 1
  
  % load tables
  if ~existfile('measurementTable.mat')
    error('Measurement table does not exist.');
  end
  load('measurementTable.mat');
  
  % loop through obsTimes
  for i = 1:length(obsTimes)
    currentTime = obsTimes(i);
    
    % loop through all tables
    W{i} = [];
    H{i} = [];
    measurementIndexTable{i} = []; %#ok<*AGROW>
    for j = 1:length(sigmaVectors) %#ok<USENS>
      entry = sigmaVectors{j};
      index = find(entry(:,1) == currentTime);
      if ~isempty(index)
        sigma = entry(index,2);
        W{i} = [W{i},sigma^2];
        vec = sparse(zeros(1,numVar));
        vec(end-numMeas+measurementIndex(j)) = 1; %#ok<*SPRIX>
        H{i} = [H{i};vec];
        measurementIndexTable{i} = [measurementIndexTable{i},...
                                    measurementIndex(j)];
      end
    end
    W{i} = diag(W{i});
    if isempty(W{i})
      warning('rfmatlab:eclipseKalman',['Something wrong with either measurement table or ' ...
               'timespec. ', 'Time ',num2str(currentTime),' not present.']);
    end
  end
  
% build H and W using predefined variables
else 
  
  % load input data
  if existfile('inputData.mat')
    load('inputData.mat');
  end
  
  % build H
  if ~exist('H','var')
    %H = sparse(zeros(numMeas,numVar+numMeas));
    H = sparse(zeros(numMeas,numVar));
    for i = 1:numMeas
      H(i,end-numMeas+i) = 1;
    end
  end
  
  % maybe not all measurements should be used
  if isfield(kalmanOptions,'Hones') && ...
             kalmanOptions.Hones==1
      c = H(:,2);
  else
      [~,c] = find(H~=0);
  end
  c = c - numVar + numMeas;
  numMeas = length(c);
  measurements = measurements(c,:);
  
  % build W
  if ~exist('W','var')
    W = eye(numMeas,numMeas);
  end
  
  % check if absolute weight is specified
  if isfield(kalmanOptions,'absolute')
    index = kalmanOptions.absolute;
    if ~isfield(kalmanOptions,'absoluteWeight')
      error('Absolute weight not specified.');
    end
    absoluteWeight = kalmanOptions.absoluteWeight;
    if length(absoluteWeight) == 1
      absoluteWeight = absoluteWeight*ones(1,length(index));
    end
    if length(index) ~= length(absoluteWeight)
      error('Wrong dimention of absolute weight.');
    end
    W(index,index) = diag(absoluteWeight);
  end
  
  % save data
  save('inputData.mat','-append','H','W');
  
  % copy matrices to cell arrays
  Htilde = H;
  Wtilde = W;
  clear H W;
  for i = 1:length(obsTimes)
    H{i} = sparse(Htilde);
    W{i} = sparse(Wtilde);
    measurementIndexTable{i} = 1:size(Htilde,1);
  end
  
  % check if relative weight is specified
  if isfield(kalmanOptions,'relative')
    index = kalmanOptions.relative;
    if ~isfield(kalmanOptions,'relativeWeight')
      error('Relative weight not specified.');
    end
    relativeWeight = kalmanOptions.relativeWeight;
    if length(relativeWeight) == 1
      relativeWeight = relativeWeight*ones(1,length(index));
    end
    if length(index) ~= length(relativeWeight)
      error('Wrong dimention of relative weight.');
    end
    for i = 1:length(obsTimes)
      W{i}(index,index) = sparse(diag(1e-3 + (abs(measurements(index,i)).*...
                                       relativeWeight).^2));
    end
  end
  
end
