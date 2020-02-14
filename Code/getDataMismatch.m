function [obj,objStd,objReal,varargout] = getDataMismatch(simData,W,measurement)

% calculate data mismatch of the ensemble 
% function [obj,objStd,objReal,varargout]=getDataMismatch(simData,W,measurment)

% INPUTS:
% ------
% simData: 			simulated measurments (normalized by the stds of observation errors) 
% W:				weight vector for observation elements (all equal to 1 with the normalization)
% measurement:		real measurments (normalized by the stds of observation errors) 		

% OUTPUTS:
% -------
% obj: 				average data mismatch (scalar)
% objStd:			std of data mismatch (scalar)
% objReal: 			data mismatch of the ensemble of models (vector)
% varargout (mismatchMtx): 	ensemble of data mismatch w.r.t each observation element (matrix)			
%
%

% Copyright (c) 2010-2014 IRIS, All Rights Reserved.

% $Id: //depot/rfmatlab/main/eclipseKalman/getDataMismatch.m#7 $

% $DateTime: 2017/05/18 15:13:28 $

ne=size(simData,2);
objReal=zeros(ne,1);
mismatchMtx = zeros(size(simData)); 

for j=1:ne
    mismatchMtx(:,j) = ((simData(:,j)-measurement).^2)./(W.^2);
    objReal(j)=sum(mismatchMtx(:,j)); % for diagonal weight matrices only
end
% mean squared error over the ensemble
obj=mean(objReal);
objStd=std(objReal);

if nargout >= 4
    varargout{1} = mismatchMtx;
end
if nargout >= 5
    load('trueSolutionSmoother.mat','num_prod','num_seis');
    if size(simData,1) ~= (num_prod + num_seis)
        error('getDataMismatch_prodSeis: Measurement sizes do not match');
    end

    prodMismatchMtx= mismatchMtx(1:num_prod,:);
    seisMismatchMtx = mismatchMtx((1+num_prod):end,:);
    for j=1:ne
        objReal_prod(j)=sum(prodMismatchMtx(:,j)); %#ok<*AGROW> % for diagonal weight matrices only
        objReal_seis(j)=sum(seisMismatchMtx(:,j)); % for diagonal weight matrices only
    end
    % mean squared error over the ensemble
    obj_prod=mean(objReal_prod);
    objStd_prod=std(objReal_prod);
    obj_seis=mean(objReal_seis);
    objStd_seis=std(objReal_seis);
    varargout{2} = {objReal_prod,objReal_seis,[obj_prod objStd_prod],[obj_seis objStd_seis]};
end
if nargout > 5
    error('getDataMismatch_prodSeis: too many output arguments');
end
