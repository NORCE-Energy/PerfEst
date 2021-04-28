function ensemble=generateInitialEnsemble(kalmanOptions,options)

% function ensemble=generateInitialEnsemble(kalmanOptions,options)
%
% Generate an initial ensemble.
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/generateInitialEnsemble.m#1 $
% $DateTime: 2011/10/13 14:36:15 $


% size of the field, used for static variables:
fieldSize=options.fieldSize;

% Either read the inital ensemble from file specified in
% kalmanOptions.initialEnsemble or generate it: 
if isfield(kalmanOptions,'initialEnsemble')
    load(kalmanOptions.initialEnsemble,'ensemble')
else
  % initialization of ensemble, that is, initialization of the
  % static variables:
  ensemble=zeros(length(kalmanOptions.staticVarMean), ...
		 kalmanOptions.ensembleSize); 
  if isfield(kalmanOptions,'state')
      s = 309+kalmanOptions.state;
      while s > 4.2897e+09, s = round(s/1000); end
      rng(s); 
  end
  if 0 %fieldSize > 1 % add condition here, given from kalmanOptions, for now it is hard coded.
      for i=1:kalmanOptions.ensembleSize
          numStart=1;
          numEnd=fieldSize;
          for k=1:size(options.staticVar,1)
              [ensemble(numStart:numEnd,i),~,pert]= ...
                  gaussianWithVariableParameters(options.dim, ...
                  kalmanOptions.staticVarMean(numStart:numEnd), ...
                  kalmanOptions.staticVarStdDev(numStart:numEnd), ...
                  kalmanOptions.meanCorrLength,...
                  kalmanOptions.stdCorrLength,kalmanOptions);
              if strcmp(deblank(options.staticVar(k,:)),'POROART')
                  pert1=pert;
              elseif strcmp(deblank(options.staticVar(k,:)),'POROVEN')
                  pert2=pert;
              end
              if strcmp(deblank(options.staticVar(k,:)),'POROVEN')
                  ensemble(numStart:numEnd,i)=kalmanOptions.staticVarMean(numStart:numEnd)+...
                      kalmanOptions.staticVarStdDev(numStart:numEnd).*(-0.8*pert1+0.6*pert2);
              end
              numStart=numStart+fieldSize;
              numEnd=numEnd+fieldSize;
          end
      end
  else
      for i=1:kalmanOptions.ensembleSize
          numStart=1;
          numEnd=fieldSize;
          for k=1:size(options.staticVar,1)
              [ensemble(numStart:numEnd,i)]= ...
                  gaussianWithVariableParameters(options.dim, ...
                  kalmanOptions.staticVarMean(numStart:numEnd), ...
                  kalmanOptions.staticVarStdDev(numStart:numEnd), ...
                  kalmanOptions.meanCorrLength,...
                  kalmanOptions.stdCorrLength,kalmanOptions);
              numStart=numStart+fieldSize;
              numEnd=numEnd+fieldSize;
          end
      end
  end
      
  
  % only keep active cells
  if isequal(size(ensemble,1),options.fieldSize* size(options.staticVar,1))
    disp('ensemble should only contain values for active cells')
    disp('this is fixed for you')
    numVar=0;
    if isfield(options,'freeparam')
        numVar=sum(options.freeparam);
    end
        fullEnsemble=ensemble;
        ensemble=zeros(options.numGridBlocks*size(options.staticVar,1),kalmanOptions.ensembleSize);
    for i=1:size(options.staticVar,1)
        ensPart= ...
            fullEnsemble((numVar+1+(options.fieldSize*(i-1))):(numVar+(options.fieldSize*i)),:);
        ensemble((numVar+1+(options.numGridBlocks*(i-1))):(numVar+(options.numGridBlocks*i)),:)=...
            ensPart(options.actnum==1,:);
    end
    if isfield(kalmanOptions,'initialEnsemble')
        save(kalmanOptions.initialEnsemble,'ensemble');
    end
  end

end

% if bounds are specified, make sure that they are respected:
numStaticVar=1:size(options.staticVar,1)*options.numGridBlocks;

if isfield(kalmanOptions,'freeparamLB')
    ensemble(1:sum(options.freeparam),:)= ...
      adjustVariableWithInBounds(ensemble(1:sum(options.freeparam),:), ...
				 kalmanOptions.freeparamLB, ...
				 kalmanOptions.freeparamUB);
    numStaticVar=numStaticVar+sum(options.freeparam);
end
% adjust staticVar within bounds if required:
if isfield(kalmanOptions,'staticVarUB')
  ensemble(numStaticVar,:)= ...
      adjustVariableWithInBounds(ensemble(numStaticVar,:), ...
				 kalmanOptions.staticVarLB, ...
				 kalmanOptions.staticVarUB);
end

end
