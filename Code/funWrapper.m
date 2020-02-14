function  [Y,volTracerArt,volTracerCap,volTracerVen,velX,velY,velZ,rateQ,pres,modifiedOptions]=funWrapper(options,param)

%setupOptions; % default options
if nargin >1
    if isfield(options,'variableList')
        index=1;
        for I=1:size(options.variableList,1)
            if options.variableLength(I)==1
                options.(strip(options.variableList(I,:)))=param(index:index+options.variableLength(I)-1)*options.(strip(options.variableList(I,:)));
            else
                options.(strip(options.variableList(I,:)))=...
                    reshape(param(index:index+options.variableLength(I)-1),size(options.(strip(options.variableList(I,:))))).*options.(strip(options.variableList(I,:)));
            end
            index=index+options.variableLength(I);
        end
    end
end

options.noFlowSouthQ1 = 1;
options.noFlowNorthQ1 = 1;
options.noFlowSouthQ2 = 1;
options.noFlowNorthQ2 = 1;

modifiedOptions=options;
[volTracerArt,volTracerCap,volTracerVen,velX,velY,velZ,rateQ,pres]=runFlowSolver3D2QSort(options);

% calculate concentrations per block
timeSize=size(volTracerArt,2);
Y=volTracerArt.*repmat(options.porosityArt(:),1,timeSize)+...
    volTracerCap.*repmat(options.porosityQ(:),1,timeSize)+...
    volTracerVen.*repmat(options.porosityVen(:),1,timeSize);
%Y(:,1)=[]; % at zero...
%Y=Y(:,10:100:end);
for I=1:length(options.reportTime)
    [~,b]=min(abs(options.reportTime(I)-options.time));
    v(I)=b;
end
Y=Y(:,v);
%Y=Y(:,(size(Y,2)-1)/60:size(Y,2));
return