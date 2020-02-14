%load problemCase
%load testCase
load modOpt
load('../classicDecon/aif')

modifiedOptions.time=0:0.1:90;
[Y,volTracerArt,volTracerCap,volTracerVen,velX,velY,velZ,rateQ,pres,modifiedOptions]=funWrapper(modifiedOptions);
plot(sum(Y)),title('total amount of contrast in system')

col=aif(1:89);
row=zeros(1,89);
row(1,1)=col(1,1);
mtxA=toeplitz(col,row);
[U,S,V]=svd(mtxA);
% figure
% hold on
% plot(V(:,1),'LineWidth',2)
% plot(V(:,5),'LineWidth',2)
% plot(V(:,10),'LineWidth',2)
% hold off

% cTot=sum(meas,1);
% coeff=cTot*U(:,:);
% diagS=diag(S);
% coeffSinv=coeff'./diagS(:);
% figure('Name','Total residual function')
%     hold on
% plot(V(:,1:38)*coeffSinv(1:38),'b','LineWidth',2)
% plot(V(:,1:50)*coeffSinv(1:50),'g','LineWidth',2)
% plot(V(:,1:52)*coeffSinv(1:52),'r','LineWidth',2)
% plot(V(:,1:64)*coeffSinv(1:64),'k','LineWidth',2)
% legend('38','50','52','64');

resMax=zeros(size(Y,1),1);
for jj=1:size(Y,1)
    if modifiedOptions.mask(jj)
        coeff=Y(jj,:)*U(:,:);
        diagS=diag(S);
        coeffSinv=coeff'./diagS(:);
        flowScaledRes=V(:,1:50)*coeffSinv(1:50);
        resMax(jj)=max(flowScaledRes);
    end
end

AvgPerf=sum(resMax(61*66+1:end))/sum(modifiedOptions.mask(61*66+1:end))

figure('Name','rateQ vs max residual, k(t)')
plot(rateQ(61*66+1:end),resMax(61*66+1:end),'.');