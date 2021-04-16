function [x,corrLength,pert]=gaussianWithVariableParameters(fieldDim,meanValue, ...
			  Sdev,meanCorrLength,stdCorrLength,options)

% function [x,corrLength,pert]=gaussianWithVariableParameters(meanValue, ...
%			  Sdev,meanCorrLength,stdCorrLength); 
%
% Setup a gaussian field with correlation length drawn from a
% normal distribution.
%
% x              : the generated field.
% corrLenght     : the drawn correlation length (optional). 
% pert           : pertubation used to generate the field (optional).
%
% fieldDim       : dimension of the field.
% meanValue      : the mean value of the field (vector of the size
%                  of the field).
% Sdev           : standard deviation of the field.
% meanCorrLength : mean correlation length.
% stdCorrLength  : standard deviation of the correlation length.
% options        :
%  -.indepHorLayers: the horizontal layers are genereated
%                    independently if this field exist. 
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/gaussianWithVariableParameters.m#1 $
% $DateTime: 2011/10/13 14:36:15 $


corrLength=addgnoise(meanCorrLength,stdCorrLength,1);
if length(fieldDim)<3
  x=meanValue+fastGaussian(fieldDim,Sdev,corrLength); 
else
  if length(corrLength)==2 || (nargin > 5 && isfield(options, ...
						   'indepHorLayers')) 
    layerDim=prod(fieldDim(1:2));
    % initialization:
    x=meanValue;
    if length(Sdev)==1
      for I=1:fieldDim(3)
	x(1+(I-1)*layerDim:I*layerDim,1)=meanValue(1+(I-1)*layerDim: ...
						 I*layerDim)+ ...
	    fastGaussian(fieldDim(1:2),Sdev,corrLength);   
	% generate new correlation length for the next layer:
	corrLength=addgnoise(meanCorrLength,stdCorrLength,1);
      end
    else
      for I=1:fieldDim(3)
	x(1+(I-1)*layerDim:I*layerDim,1)=meanValue(1+(I-1)*layerDim: ...
						 I*layerDim)+ ...
	    fastGaussian(fieldDim(1:2),Sdev(1+(I-1)*layerDim:I* ...
				       layerDim),corrLength);    
	% generate new correlation length for the next layer:
	corrLength=addgnoise(meanCorrLength,stdCorrLength,1);
      end
    end
  else
      if nargout>2
          pert=fastGaussian3d(fieldDim,Sdev,corrLength);
          x=meanValue+pert;
          pert=pert./Sdev;
      else
          x=meanValue+fastGaussian3d(fieldDim,Sdev,corrLength);
      end
  end
end


