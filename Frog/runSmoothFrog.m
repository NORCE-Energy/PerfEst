function runSmoothFrog

if ~exist('initialState.mat','file')
    initializeFrog;
end
load initialState.mat;

% Iterative ensemble smoother
runIES;

% Load final ensemble
fileName = @(x) ['ensemble',num2str(x),'.mat'];
A = dir(fileName('*'));
iter = length(A)-1;
if iter <= 0, error('Nothing has happend'); end
load(['ensemble',num2str(iter),'.mat']);
xValue = mean(ensemble,2);
options = setOptions(options,xValue); %#ok<*NODEF>

%Post-process ...
state=runFlowSolver3D2QSort(options); %#ok<*NASGU>
save('finalState.mat','state','options');

% Plot results
plotResults;


   
    
    
    
   
    



  
  
  
  
  
