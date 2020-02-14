function [mm,cmi,mmit,lambdaiterhist,lambdahist,it,totalit,objred,InvErr,obj,objd,objm] = LM_1c(m_pr,cmi,d_obs,cdi,options,dfunc,Jfunc,m_pr0,varargin)
% Levenberg-Marquardt algorithm for finding a  minimum
% Used for  seismic inversion

% Input:
% m_pr: prior mean
% cmi: inverse of prior covariance
% d_obs: observed  data
% cdi: inverse of measurement covariance. Column vector or number
% options
% dfunc: name of function that calculates the data
% Jfunc: name of function that calculates the sensitivity  matrix
% m_pr0: original prior (only used in stop criteria for doublediff.)
% varargin: variables needed by dfunc and  Jfunc

% Output:
% mm: MAP estimate
% cmi: inverse of posterior covariance
% mmit: all accepted steps
% lambdaiterhist: vector with number of lambdaiterations
% lambdahist: vector with lambda
% it:  number of outer iterations
% totalit: number of iterations
% objred: reduction in obj
% InvErr: error compared to true
% obj: total objective function
% objd: data objective function
% objm: model objective function


mm=m_pr;
N=length(mm);
lambdaWasIncreased=0;
totalit=0;
min_InvErr=1;
mmit=[];
if isfield(options,'IT_max')
    IT_max = options.IT_max;
else
    IT_max =20;
end

if isfield(options,'maxLambdaIt')
    maxLambdaIt = options.maxLambdaIt;
else
    maxLambdaIt =5;
end

if isfield(options,'lambda')
    lambda = options.lambda;
else
    lambda =100;
end

if isfield(options,'lambdafactor')
    lambdafactor = options.lambdafactor;
else
    lambdafactor =10;
end

if isfield(options,'lambdaMin')
    lambdaMin = options.lambdaMin;
else
    lambdaMin =0.01;
end

if isfield(options,'doubledifference')
    doubledifference = options.doubledifference;
else
    doubledifference =0;
end

if isfield(options,'minReduction')
    minReduction = options.minReduction;
else
    minReduction =0.2;
end

if isfield(options,'dxmin')
    dxmin = options.dxmin;
else
    dxmin =0.01;
end

if isfield(options,'trueSolution') 
    synthetic = 1;
    true=options.trueSolution;
else
    synthetic =0;
end


%calculate data and obj func
[d_calc,T] =feval( dfunc,mm,varargin );

objd=getMismatch(d_calc,cdi,d_obs);
objm=0;
obj=objd+objm;

disp(['obj=',num2str(obj)]);
disp(['objd=',num2str(objd)]);
disp(['objm=',num2str(objm)]);


it=0;
while it < IT_max
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    it = it + 1;
        
    %calculate sensitivity
    Jb=feval(Jfunc,m_pr,varargin,T);
    
    %observed minus simulated data
    d_dr = d_obs -d_calc;
   
    if it==1
        cmtiOld=cmi;
    else
        cmtiOld=cmti;
    end
    
    if size(cdi,1)>1
        cmti=cmi+real(bsxfun(@times,Jb',cdi')*Jb);
    else
        cmti=cmi+cdi(1)*real((Jb'*Jb));
    end
    
    iterLambda=1;
    smallReduction=0;
    while iterLambda<maxLambdaIt
        mmold=mm;
        dataOld=d_calc;
        
        %update
        mm = mm + (cmti+lambda*eye(N))\real(Jb'*(cdi.*d_dr)+cmi*(m_pr-mm));
        
        %calculate data and obj func
        [d_calc,T] =feval( dfunc,mm,varargin );
        
        objdnew=getMismatch(d_calc,cdi,d_obs);
        objmnew=getMismatch(mm,cmi,m_pr);
        objnew=objdnew+objmnew;
        
        disp(['obj=',num2str(obj)]);
        disp(['objnew=',num2str(objnew)]);
        disp(['objdnew=',num2str(objdnew)]);
        disp(['objmnew=',num2str(objmnew)]);
        
        %if objective function has decreased, the step is accepted.
        %If not, the calculation is redone with a larger lambda.
        if objnew>obj
            lambda=lambda*lambdafactor;
            if iterLambda>1
                lambda=lambda*lambdafactor;%increase more if two iterations have already been done
            end
            disp(['increasing Lambda to ',num2str(lambda)]);
            iterLambda=iterLambda+1;
            mm=mmold;
            d_calc=dataOld;
            lambdaWasIncreased=1;
        else
            lambdahist(it)=lambda;
            if  lambda>lambdaMin
                if lambdaWasIncreased==0
                    lambda=lambda/lambdafactor;
                    disp(['reducing Lambda to ',num2str(lambda)]);
                else
                    lambdaWasIncreased=0;
                    disp(['Lambda is kept at ',num2str(lambda)]);
                end
            else
                disp(['Lambda is kept at ',num2str(lambda)]);
                lambdaWasIncreased=0;
            end
            reduc=abs(objnew-obj)/abs(obj)*100;
            disp(['reduction in obj for this iter. in percent: ',num2str(reduc)])
            if reduc<minReduction
                smallReduction=1;
            end
            objred(it)=reduc;
            
            obj=objnew;
            objd=objdnew;
            objm=objmnew;
            mmit=[mmit;mm];
            break % break the inner loop over lambda
        end
        
    end
    totalit=totalit+iterLambda;
    disp(['Number of iterations=',num2str(totalit)]);
    lambdaiterhist(it)=iterLambda;
    
 
    %change in model variable
    if doubledifference==0
        dx=norm(mmold  - mm)/norm(mmold)
    else
        dxNotDoubleDiff=norm(mmold  - mm)/norm(mmold)
        dx=norm(mmold  - mm)/(norm(mmold-m_pr0)+1)
    end
    
    mm_mat(1:N, it) = mm;
    
    %compare with true if synthetic
    if synthetic
        InvErr(it)=norm(mm-true)/norm(true);
        
        if InvErr(it) < min_InvErr
            min_InvErr = InvErr(it);
        end
        InvError = InvErr(it)
        MinimumInvError = min_InvErr
    end
    
    %check stopping criteria
    if smallReduction==1 && dx<dxmin && it>1
        cmi=cmti;
        tline = ['terminating iterations: reduction of objective function is less than ', ...
            num2str(minReduction),'%'];
        disp(tline);
        break
    end
    if iterLambda==maxLambdaIt
        cmi=cmtiOld;
        tline = 'terminating iterations: iterLambda>=maxInnerIter';
        disp(tline);
        break
    end
    if it==IT_max
        tline = 'terminating iterations: it>=IT_max';
        disp(tline);
        cmi=cmti;
    end
end



end
