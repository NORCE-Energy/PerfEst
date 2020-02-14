
% create options for runFlowSolver3D2Q
clear options

%prm=rmfield(prm,'etafunc');
% here is the mask:

% time
options.time=0:0.1:6;
%options.time=0:0.1:60;
%options.time=0:0.001:60;
%Discretization
nx = 10;
ny = 10;
nz = 10;

%Domain
xL = 0.01;
yL = 0.01;
zL = 0.01;

%Viscosity
visc = 0.004*ones(nx,ny,nz,2); % 4*1e-3*ones(nx,ny,nz,2); % the value is 3e-6 here..

%Permeability in x-direction
permX = 2e-10*ones(nx,ny,nz,2);
%Permeability in y-direction
permY = 2e-10*ones(nx,ny,nz,2);
%Permeability in z-direction
permZ = 2e-10*ones(nx,ny,nz,2);

permX(:,:,:,2)=2e-10;
permY(:,:,:,2)=2e-10;
permZ(:,:,:,2)=2e-10;

%Inter-compartment conductivities
%This is a dimensionless quantity: Permeability [m^2] divided by
%interfacial area in a one one unit volume
%permQ = 2e-10*ones(nx,ny,nz);
permQ = 1e-8*ones(nx,ny,nz);
% porosities
porosityArt=0.05*ones(nx,ny,nz);
porosityVen=0.11*ones(nx,ny,nz);
porosityQ=0.02*ones(nx,ny,nz);


% inflow profile
cValBnd = 5;%25e6;
timeStart=0;
cValBndVen=0;%25e6;
timeStartVen=6;


% boundary conditions
prsValueWestQ1=133*85;
prsValueEastQ1=133*85;
prsValueSouthQ1=133*85;
prsValueNorthQ1=133*85;
prsValueDownQ1=133*85;
prsValueUpQ1=133*85;
prsValueWestQ2=133*10;
prsValueEastQ2=133*10;
prsValueSouthQ2=133*10;
prsValueNorthQ2=133*10;
prsValueDownQ2=133*10;
prsValueUpQ2=133*10;

% no-flow up and down - note that this overides any corresponding prs cnd
options.noFlowWestQ1 = 1;
options.noFlowEastQ1 = 1;
options.noFlowSouthQ1 = 1;
options.noFlowNorthQ1 = 1;
options.noFlowDownQ1 = 0;
options.noFlowUpQ1 = 1;
options.noFlowDownQ2 = 0;
options.noFlowUpQ2 = 1;
options.noFlowEastQ2=1;
options.noFlowWestQ2=1;
options.noFlowSouthQ2 = 1;
options.noFlowNorthQ2 = 1;


options.xL=xL;
options.yL=yL;
options.zL=zL;
options.nx=nx;
options.ny=ny;
options.nz=nz;
options.visc=visc;
options.permX=permX;
options.permY=permY;
options.permZ=permZ;
options.permQ=permQ;
options.porosityArt=porosityArt;
options.porosityQ=porosityQ;
options.porosityVen=porosityVen;

options.prsValueWestQ1 = prsValueWestQ1;
options.prsValueEastQ1 = prsValueEastQ1;
options.prsValueSouthQ1 = prsValueSouthQ1;
options.prsValueNorthQ1 = prsValueNorthQ1;
options.prsValueDownQ1 = prsValueDownQ1 ;
options.prsValueUpQ1 = prsValueUpQ1;
options.prsValueWestQ2 = prsValueWestQ2;
options.prsValueEastQ2 = prsValueEastQ2;
options.prsValueSouthQ2 = prsValueSouthQ2;
options.prsValueNorthQ2 = prsValueNorthQ2;
options.prsValueDownQ2 = prsValueDownQ2;
options.prsValueUpQ2 = prsValueUpQ2;
options.porosityArt=porosityArt;
options.porosityVen=porosityVen;
options.porosityQ=porosityQ;
options.cValBnd=cValBnd;
options.timeStart=timeStart;
options.gammaShape=1.1;
options.gammaScale=1.5;
%options.timeEnd=timeEnd;

options.mask = ones(nx,ny,nz,2);
options.maskQ = ones(nx,ny,nz);
options.mask(:,:,1,:) = 0;
options.mask(5:6,5:6,1,1) = 1;
options.mask(1:9:10,1:9:10,1,2) = 1;

%Perm pertubation: (base case)

% % % options.permX(5:6,5:6,1:7,1) = 0.1*options.permX(5:6,5:6,1:7,1);
% % % options.permY(5:6,5:6,1:7,1) = 0.1*options.permY(5:6,5:6,1:7,1);
% % % options.permZ(5:6,5:6,1:7,1) = 10*options.permZ(5:6,5:6,1:7,1);
% % % 
% % % options.permX(3:8,5:6,8,1) = 10*options.permX(3:8,5:6,8,1);
% % % options.permY(3:8,5:6,8,1) = 10*options.permY(3:8,5:6,8,1);
% % % options.permZ(3:8,5:6,8,1) = 10*options.permZ(3:8,5:6,8,1);
% % % 
% % % options.permX(5:6,7:8,8,1) = 10*options.permX(5:6,7:8,8,1);
% % % options.permY(5:6,7:8,8,1) = 10*options.permY(5:6,7:8,8,1);
% % % options.permZ(5:6,7:8,8,1) = 10*options.permZ(5:6,7:8,8,1);
% % % 
% % % options.permX(5:6,3:4,8,1) = 10*options.permX(5:6,3:4,8,1);
% % % options.permY(5:6,3:4,8,1) = 10*options.permY(5:6,3:4,8,1);
% % % options.permZ(5:6,3:4,8,1) = 10*options.permZ(5:6,3:4,8,1);

options.maskQ(5:6,5:6,1:7) = 0;
options.maskQ(3:8,5:6,8) = 0;
options.maskQ(5:6,7:8,8) = 0;
options.maskQ(5:6,3:4,8) = 0;


  
load('tofArt.mat','tofArtObserved','permXBase','permYBase','permZBase');

%Teste baseløsning ...
% options.permX=permXBase;
% options.permY=permYBase;
% options.permZ=permZBase;

[state]=runFlowSolver3D2QSort(options);

%Lagre baseløsning ...
% % % % % %   tofArtObserved=state.tofArt;
% % % % % %   permXBase = options.permX;
% % % % % %   permYBase = options.permY;
% % % % % %   permZBase = options.permZ;
% % % % % %   save('tofArt.mat','tofArtObserved','permXBase','permYBase','permZBase');

itMax=12;
d_perm_log = zeros(itMax+1,1);
grad_perm_log = zeros(itMax+1,1);
tot_res_log = zeros(itMax+1,1);
delta_log = zeros(itMax+1,1);
d_t_log = zeros(itMax+1,1);
deltaStart=-4.0;
delta_reduction_factor = 0.5;

alpha=0.0;
permRef=2e-10;
alphaPermReg = alpha/permRef;

  
figure('Name','kx/ky/kz','NumberTitle','off')
hold on
colormap(jet(128));

for itCnt=1:itMax
  [gradPerm, gradPor, sampleMask, lamTau, lamSrc, lamPrs] = computGradAdj(options, state, tofArtObserved, alphaPermReg, 0);
  
  if itCnt==1
    d_tofArt=(state.tofArt-tofArtObserved).*reshape(sampleMask,nx,ny,nz);
    d_t = 0.5*reshape(d_tofArt,1,[])*reshape(d_tofArt,1,[])';
    d_t_log(1)=d_t;

    d_permX = norm(reshape(options.permX-permXBase,1,[])); % .*sampleMask);
    d_permY = norm(reshape(options.permY-permYBase,1,[])); % .*sampleMask);
    d_permZ = norm(reshape(options.permZ-permZBase,1,[])); % .*sampleMask);
    d_perm = sqrt(d_permX*d_permX+d_permY*d_permY+d_permZ*d_permZ);
    d_perm_log(1)=d_perm;
    
    totalResidReg = d_t;
    
    %regularization terms
%     totalResidReg = totalResidReg + (alphaPermReg*norm(reshape(options.permX,1,[]).*sampleMask))^2;
%     totalResidReg = totalResidReg + (alphaPermReg*norm(reshape(options.permY,1,[]).*sampleMask))^2;
%     totalResidReg = totalResidReg + (alphaPermReg*norm(reshape(options.permZ,1,[]).*sampleMask))^2;
    tot_res_log(1)=totalResidReg;
    
    for j=1:10
      climsTofErr = [-0.2 0.2];
      subplot(10,10,10*(j-1)+7);
      imagesc(d_tofArt(:,:,j),climsTofErr);
      if j==1
        c=colorbar;
        c.Location='north';
      end
      ax = gca;
      ax.XAxis.Visible = 'off';
      ax.YAxis.Visible = 'off';

      climsLamTau = [-1.0e9 1.0e9];
      subplot(10,10,10*(j-1)+8);
      imagesc(lamTau(:,:,j),climsLamTau);
      if j==1
        c=colorbar;
        c.Location='northoutside';
      end
      ax = gca;
      ax.XAxis.Visible = 'off';
      ax.YAxis.Visible = 'off';   

      climsLamSrc = [-30 30];
      subplot(10,10,10*(j-1)+9);
      imagesc(lamSrc(:,:,j),climsLamSrc);
      if j==1
        c=colorbar;
        c.Location='northoutside';
      end
      ax = gca;
      ax.XAxis.Visible = 'off';
      ax.YAxis.Visible = 'off';   

      climsLamPrs = [-10e11 -5e11];
      subplot(10,10,10*(j-1)+10);
      imagesc(lamPrs(:,:,j,1),climsLamPrs);
      if j==1
        c=colorbar;
        c.Location='northoutside';
      end
      ax = gca;
      ax.XAxis.Visible = 'off';
      ax.YAxis.Visible = 'off';   


      set(gca,'YDir','normal')
    end
    
    drawnow,pause(1)
    
  end
  
  %regularization terms - now included in computGradAdj ...
%   gradPerm(:,:,:,1,1) = gradPerm(:,:,:,1,1) + alpha*(options.permX(:,:,:,1) - permRef);
%   gradPerm(:,:,:,1,2) = gradPerm(:,:,:,1,2) + alpha*(options.permY(:,:,:,1) - permRef);
%   gradPerm(:,:,:,1,3) = gradPerm(:,:,:,1,3) + alpha*(options.permZ(:,:,:,1) - permRef);

  %sampleMask??? - now included in computGradAdj ...
% %   gradPermX = reshape(reshape(gradPerm(:,:,:,1,1),1,[]).*sampleMask',nx,ny,nz);
% %   gradPermY = reshape(reshape(gradPerm(:,:,:,1,2),1,[]).*sampleMask',nx,ny,nz);
% %   gradPermZ = reshape(reshape(gradPerm(:,:,:,1,3),1,[]).*sampleMask',nx,ny,nz);
  
  gradPermX = reshape(reshape(gradPerm(:,:,:,1,1),1,[]),nx,ny,nz);
  gradPermY = reshape(reshape(gradPerm(:,:,:,1,2),1,[]),nx,ny,nz);
  gradPermZ = reshape(reshape(gradPerm(:,:,:,1,3),1,[]),nx,ny,nz);
  
  permX_old = options.permX(:,:,:,1);
  permY_old = options.permY(:,:,:,1);
  permZ_old = options.permZ(:,:,:,1);
  
  % log perm ...
  gradLogPermX = gradPermX .* permX_old;
  gradLogPermY = gradPermY .* permY_old;
  gradLogPermZ = gradPermZ .* permZ_old;
  
  gradLogPermXAbsSq = reshape(gradLogPermX,1,[])*reshape(gradLogPermX,1,[])';
  gradLogPermYAbsSq = reshape(gradLogPermY,1,[])*reshape(gradLogPermY,1,[])';
  gradLogPermZAbsSq = reshape(gradLogPermZ,1,[])*reshape(gradLogPermZ,1,[])';

  gradLogPermAbsSq=gradLogPermXAbsSq+gradLogPermYAbsSq+gradLogPermZAbsSq+1.0e-25;
  
  gX=sqrt(gradLogPermXAbsSq/gradLogPermAbsSq);
  gY=sqrt(gradLogPermYAbsSq/gradLogPermAbsSq);
  gZ=sqrt(gradLogPermZAbsSq/gradLogPermAbsSq);
  
  gradLogPermX = (gX>0.0) * gradLogPermX / sqrt(gradLogPermAbsSq);
  gradLogPermY = (gY>0.0) * gradLogPermY / sqrt(gradLogPermAbsSq);
  gradLogPermZ = (gZ>0.0) * gradLogPermZ / sqrt(gradLogPermAbsSq);
  
  dfdk1Art=gradLogPermX
  dfdk2Art=gradLogPermY
  dfdk3Art=gradLogPermZ
  % ... log perm
  
  % no-log perm ...
% % %   gradPermXAbsSq = reshape(gradPermX,1,[])*reshape(gradPermX,1,[])';
% % %   gradPermYAbsSq = reshape(gradPermY,1,[])*reshape(gradPermY,1,[])';
% % %   gradPermZAbsSq = reshape(gradPermZ,1,[])*reshape(gradPermZ,1,[])';
% % % 
% % %   gradPermAbsSq=gradPermXAbsSq+gradPermYAbsSq+gradPermZAbsSq+1.0e-25;
% % %   
% % %   gradPermX = gradPermX / sqrt(gradPermAbsSq);
% % %   gradPermY = gradPermY / sqrt(gradPermAbsSq);
% % %   gradPermZ = gradPermZ / sqrt(gradPermAbsSq);
  % ... no-log perm
  
  delta = deltaStart;
  
  msg=[' gradLogPermX: ', num2str(gX),' gradLogPermY: ',num2str(gY),' gradLogPermZ: ',num2str(gZ),' gradLogPermAbs: ',num2str(sqrt(gradLogPermAbsSq))];
  disp(msg);
  
  msg=[' itCnt: ', int2str(itCnt),' t=-cm: ',num2str(0.5*sqrt(gradLogPermAbsSq)),' f(x): ',num2str(tot_res_log(itCnt))];
  disp(msg);

  % Basic line search ...
% %   for it=1:10
% %     step=delta*permRef;
% %     stepLog = delta; %*abs(log(permRef));
% % 
% %     delta_log(itCnt+1)=delta;
% % 
% %     % log perm
% %     options.permX(:,:,:,1) = exp(log(permX_old) + stepLog * gradLogPermX);
% %     options.permY(:,:,:,1) = exp(log(permY_old) + stepLog * gradLogPermY);
% %     options.permZ(:,:,:,1) = exp(log(permZ_old) + stepLog * gradLogPermZ);
% %     % ... log perm
% % 
% %     % no-log perm ...
% % % % %       options.permX(:,:,:,1) = permX_old + step * gradPermX;
% % % % %       options.permY(:,:,:,1) = permY_old + step * gradPermY;
% % % % %       options.permZ(:,:,:,1) = permZ_old + step * gradPermZ;
% % % % % 
% % % % %       options.permX(:,:,:,1) = max(options.permX(:,:,:,1),0.01*permRef);
% % % % %       options.permY(:,:,:,1) = max(options.permY(:,:,:,1),0.01*permRef);
% % % % %       options.permZ(:,:,:,1) = max(options.permZ(:,:,:,1),0.01*permRef);
% %     % ... no-log perm
% % 
% %     [state]=runFlowSolver3D2QSort(options);
% % 
% %     d_tofArt=(state.tofArt-tofArtObserved).*reshape(sampleMask,nx,ny,nz);
% %     d_t = 0.5*reshape(d_tofArt,1,[])*reshape(d_tofArt,1,[])';
% %     totalResidReg = d_t;
% %     
% %     [gradPerm, gradPor, ~, ~, ~, ~] = computGradAdj(options, state, tofArtObserved, alphaPermReg, 0);
% %     gradPermTimesGradPermNew = reshape(gradPerm(:,:,:,1,1).*options.permX(:,:,:,1),1,[]) * reshape(gradLogPermX,1,[])' ...
% %                              + reshape(gradPerm(:,:,:,1,2).*options.permY(:,:,:,1),1,[]) * reshape(gradLogPermY,1,[])' ...
% %                              + reshape(gradPerm(:,:,:,1,3).*options.permZ(:,:,:,1),1,[]) * reshape(gradLogPermZ,1,[])';
% % 
% %     msg=[int2str(itCnt),' : ', int2str(it),' a: ',num2str(delta),' f(x+a*p): ',num2str(d_t), ...
% %                         ' : ',num2str(d_t/tot_res_log(itCnt)),' atNew: ',num2str(gradPermTimesGradPermNew)];
% %     disp(msg);
% % 
% %     %regularization terms
% % %       totalResidReg = totalResidReg + (alphaPermReg*norm(reshape(options.permX(:,:,:,1),1,[]).*sampleMask))^2;
% % %       totalResidReg = totalResidReg + (alphaPermReg*norm(reshape(options.permY(:,:,:,1),1,[]).*sampleMask))^2;
% % %       totalResidReg = totalResidReg + (alphaPermReg*norm(reshape(options.permZ(:,:,:,1),1,[]).*sampleMask))^2;
% % 
% %     if totalResidReg - 0.0*0.5*sqrt(gradLogPermAbsSq) > tot_res_log(itCnt)
% %       delta = delta*delta_reduction_factor;
% %     else
% %       break;
% %     end
% %     
% %   end %basic line search
    
  %alternative weak wolfe condition ...
  c1=0.0001;
  c2=0.9;
  alpha=0;
  beta=1e10;
  t=-deltaStart;
  for it=1:10
    
    stepLog = -t; %*abs(log(permRef));
    delta_log(itCnt+1)=-t;

    % log perm
    options.permX(:,:,:,1) = exp(log(permX_old) + stepLog * gradLogPermX);
    options.permY(:,:,:,1) = exp(log(permY_old) + stepLog * gradLogPermY);
    options.permZ(:,:,:,1) = exp(log(permZ_old) + stepLog * gradLogPermZ);
    % ... log perm

    % no-log perm ...
% % %       options.permX(:,:,:,1) = permX_old + step * gradPermX;
% % %       options.permY(:,:,:,1) = permY_old + step * gradPermY;
% % %       options.permZ(:,:,:,1) = permZ_old + step * gradPermZ;
% % % 
% % %       options.permX(:,:,:,1) = max(options.permX(:,:,:,1),0.01*permRef);
% % %       options.permY(:,:,:,1) = max(options.permY(:,:,:,1),0.01*permRef);
% % %       options.permZ(:,:,:,1) = max(options.permZ(:,:,:,1),0.01*permRef);
    % ... no-log perm

    [state]=runFlowSolver3D2QSort(options);

    d_tofArt=(state.tofArt-tofArtObserved).*reshape(sampleMask,nx,ny,nz);
    d_t = 0.5*reshape(d_tofArt,1,[])*reshape(d_tofArt,1,[])';
    totalResidReg = d_t;
    
    [gradPerm, gradPor, ~, ~, ~, ~] = computGradAdj(options, state, tofArtObserved, alphaPermReg, 0);
    gradPermTimesGradPermNew = reshape(gradPerm(:,:,:,1,1).*options.permX(:,:,:,1),1,[]) * reshape(gradLogPermX,1,[])' ...
                             + reshape(gradPerm(:,:,:,1,2).*options.permY(:,:,:,1),1,[]) * reshape(gradLogPermY,1,[])' ...
                             + reshape(gradPerm(:,:,:,1,3).*options.permZ(:,:,:,1),1,[]) * reshape(gradLogPermZ,1,[])';

    msg=[int2str(itCnt),' : ', int2str(it),' a: ',num2str(-t),' f(x+a*p): ',num2str(d_t), ...
                        ' : ',num2str(d_t/tot_res_log(itCnt)),' atNew: ',num2str(gradPermTimesGradPermNew)];
    disp(msg);
    
    if totalResidReg > tot_res_log(itCnt) - c1*t*gradLogPermAbsSq
      beta=t;
      t=0.5*(alpha+beta);
    elseif -gradPermTimesGradPermNew < -c2*gradLogPermAbsSq
      alpha=t;
      if beta==1e10
        t=2*alpha;
      else
        t=0.5*(alpha+beta);
      end
    else
      break;
    end
    
  end %weak wolfe
  
  tot_res_log(itCnt+1)=totalResidReg;
  
  %grad_perm_log(itCnt+1)=sqrt(gradPermAbsSq);
  
  delta_log(1)=0.0;

%   d_permX = norm(reshape(options.permX-permXBase,1,[])); % .*sampleMask);
%   d_permY = norm(reshape(options.permY-permYBase,1,[])); % .*sampleMask);
%   d_permZ = norm(reshape(options.permZ-permZBase,1,[])); % .*sampleMask);
%   d_perm = sqrt(d_permX*d_permX+d_permY*d_permY+d_permZ*d_permZ);
%   d_perm_log(itCnt+1)=d_perm;
  
  d_t_log(itCnt+1)=d_t
  %resid_tof_art = d_tofArt(:,:,5:5)
  
%   if mod(itCnt,1) == 0
%    semilogy(reshape(options.permX(:,5,5,1),1,[]),'DisplayName',int2str(itCnt));
%   end
  
    
  for j=1:10
    climsPerm = [-11 -5];
    subplot(10,10,10*(j-1)+1);
    imagesc(log10(options.permX(:,:,j,1)),climsPerm);
    if j==1
      c=colorbar;
      c.Location='northoutside';
    end
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    subplot(10,10,10*(j-1)+2);
    imagesc(log10(options.permY(:,:,j,1)),climsPerm);
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    subplot(10,10,10*(j-1)+3);
    imagesc(log10(options.permZ(:,:,j,1)),climsPerm);
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    
    climsGrad = [-0.5 0.5];
    subplot(10,10,10*(j-1)+4);
    imagesc(gradLogPermX(:,:,j),climsGrad);
    if j==1
      c=colorbar;
      c.Location='northoutside';
    end
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    subplot(10,10,10*(j-1)+5);
    imagesc(gradLogPermY(:,:,j),climsGrad);
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    subplot(10,10,10*(j-1)+6);
    imagesc(gradLogPermZ(:,:,j),climsGrad);
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';

    climsTofErr = [-0.2 0.2];
    subplot(10,10,10*(j-1)+7);
    imagesc(d_tofArt(:,:,j),climsTofErr);
    if j==1
      c=colorbar;
      c.Location='north';
    end
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    
    climsLamTau = [-1.0e9 1.0e9];
    subplot(10,10,10*(j-1)+8);
    imagesc(lamTau(:,:,j),climsLamTau);
    if j==1
      c=colorbar;
      c.Location='northoutside';
    end
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';   

    climsLamSrc = [-30 30];
    subplot(10,10,10*(j-1)+9);
    imagesc(lamSrc(:,:,j),climsLamSrc);
    if j==1
      c=colorbar;
      c.Location='northoutside';
    end
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';      
   
    climsLamPrs = [1e13 3e13];
    subplot(10,10,10*(j-1)+10);
    imagesc(lamPrs(:,:,j,1),climsLamPrs);
    if j==1
      c=colorbar;
      c.Location='northoutside';
    end
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';       
    
    set(gca,'YDir','normal')
    %colormap(hot)
  end
  
  drawnow,pause(1)
      
end

% semilogy(reshape(options.permX(:,5,5,1),1,[]),'LineWidth',2);
% legend('show');
% hold off

figure('Name','perm X Y Z','NumberTitle','off')
semilogy(reshape(options.permX(5,5,:,1),1,[]),'DisplayName','permX');
hold on
semilogy(reshape(options.permY(5,5,:,1),1,[]),'DisplayName','permY');
semilogy(reshape(options.permZ(5,5,:,1),1,[]),'DisplayName','permZ');
legend('show');
hold off

% figure('Name','tot_res_log','NumberTitle','off')
% semilogy(tot_res_log);

figure('Name','d_t','NumberTitle','off')
semilogy(d_t_log);

permX = options.permX(:,:,1:10,1)
permY = options.permY(:,:,1:10,1)
permZ = options.permZ(:,:,1:10,1)


  
  
  
  
  
