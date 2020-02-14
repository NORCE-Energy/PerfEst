
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

nn=nx*ny*nz;
optional.x0(1:nn) = log(options.permX(:,:,:,1));
optional.x0(nn+1:2*nn) = log(options.permY(:,:,:,1));
optional.x0(2*nn+1:3*nn) = log(options.permZ(:,:,:,1));
optional.x0=optional.x0';

optional.quitLSfail = 0;

alpha=0.0;
permRef=2e-10;
alphaPermReg = alpha/permRef;

pars.options=options;
pars.tofArtObserved=tofArtObserved;
pars.alphaPermReg=alphaPermReg;

%naive LS
itMax=25;
deltaStart=-2.0;
delta_reduction_factor = 0.5;
x0=optional.x0; %log-perm
[f0,g0] = funcGrad(x0,pars);
for itCnt=1:itMax
  p=g0/sqrt(g0'*g0);
  delta = deltaStart;
  for it=1:10
    stepLog = delta;
    x=x0+stepLog*p;
    [f,g] = funcGrad(x,pars);
    msg=[int2str(itCnt),' : ', int2str(it),' a: ',num2str(delta),' f(x+a*p): ',num2str(f), ...
                        ' : ',num2str(f/f0),' atNew: ',num2str(p'*g)];
    %disp(msg);
    if f>f0
      delta = delta*delta_reduction_factor;
    else
      x0=x;
      f0=f;
      g0=g;
      break;
    end
  end
  msg=['--- naive LS:', int2str(itCnt), ' f(x)=', num2str(f)];
  disp(msg);
end
% options.permX(:,:,:,1)=reshape(exp(x(1:nn)),nx,ny,nz);
% options.permY(:,:,:,1)=reshape(exp(x(nn+1:2*nn)),nx,ny,nz);
% options.permZ(:,:,:,1)=reshape(exp(x(2*nn+1:3*nn)),nx,ny,nz);

%bfgs
pars.fgname='funcGrad';
pars.nvar=3*nn;
optional.maxit = 25;
optional.x0=x;
%optional.nvec=10; %l-bfgs
[xValue, fValue, d] = bfgs(pars,optional);
options.permX(:,:,:,1)=reshape(exp(xValue(1:nn)),nx,ny,nz);
options.permY(:,:,:,1)=reshape(exp(xValue(nn+1:2*nn)),nx,ny,nz);
options.permZ(:,:,:,1)=reshape(exp(xValue(2*nn+1:3*nn)),nx,ny,nz);

%Post-process ...
[state]=runFlowSolver3D2QSort(options);
[gradPerm, gradPor, sampleMask, lamTau, lamSrc, lamPrs] = computGradAdj(options, state, tofArtObserved, alphaPermReg, 0);

d_tofArt=(state.tofArt-tofArtObserved).*reshape(sampleMask,nx,ny,nz)
d_t = 0.5*reshape(d_tofArt,1,[])*reshape(d_tofArt,1,[])';

d_permX = norm(reshape(options.permX-permXBase,1,[])); % .*sampleMask);
d_permY = norm(reshape(options.permY-permYBase,1,[])); % .*sampleMask);
d_permZ = norm(reshape(options.permZ-permZBase,1,[])); % .*sampleMask);
d_perm = sqrt(d_permX*d_permX+d_permY*d_permY+d_permZ*d_permZ);
  
gradPermX = reshape(reshape(gradPerm(:,:,:,1,1),1,[]),nx,ny,nz);
gradPermY = reshape(reshape(gradPerm(:,:,:,1,2),1,[]),nx,ny,nz);
gradPermZ = reshape(reshape(gradPerm(:,:,:,1,3),1,[]),nx,ny,nz);
  
% log perm ...
gradLogPermX = gradPermX .* options.permX(:,:,:,1);
gradLogPermY = gradPermY .* options.permY(:,:,:,1);
gradLogPermZ = gradPermZ .* options.permZ(:,:,:,1);
  
dfdk1Art=gradLogPermX
dfdk2Art=gradLogPermY
dfdk3Art=gradLogPermZ
    
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
  imagesc(lamSrc(:,:,j,1),climsLamSrc);
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


figure('Name','perm X Y Z','NumberTitle','off')
semilogy(reshape(options.permX(5,5,:,1),1,[]),'DisplayName','permX');
hold on
semilogy(reshape(options.permY(5,5,:,1),1,[]),'DisplayName','permY');
semilogy(reshape(options.permZ(5,5,:,1),1,[]),'DisplayName','permZ');
legend('show');
hold off

permX = options.permX(:,:,1:10,1)
permY = options.permY(:,:,1:10,1)
permZ = options.permZ(:,:,1:10,1)


  
  
  
  
  
