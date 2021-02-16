
% create options for runFlowSolver3D2Q
clear options

% 1/0 : include/skip cac calculations...
options.doCAC=1;


% time
%options.time=0:0.01:12;
options.time=0:0.05:60;
%options.time=0:0.1:60;
%options.time=0:0.001:60;
%Discretization
nx = 10;
ny = 10;
nz = 10;

%Domain
xL = 0.1;
yL = 0.1;
zL = 0.1;

% %Terminal sinks (venular):
% % % termSink.x=[0.027 0.053 0.089];
% % % termSink.y=[0.027 0.053 0.089];
% % % termSink.z=[0.027 0.053 0.089];
% % % termSink.p0=133*10*[1 1 1];
% % % termSink.perm=1e-6*[1 1 1];
%%% min test:
% termSink.x=[0.015 0.085];
% termSink.y=[0.085 0.015];
% termSink.z=[0.015 0.015];
% termSink.p0=133*10*[1 1];
% termSink.perm=1e-6*[1 1];

%%% min test:
% termSrc.x=[0.015 0.085];
% termSrc.y=[0.015 0.085];
% termSrc.z=[0.015 0.015];
% termSrc.p0=133*85*[1 1];
% termSrc.perm=1e-6*[1 1];

% BigBrain Vener: dim(346 x 448 x 319)
% % %    211   173    75
% % %    189   184   126
% % %    330   192   133
% % %    176   218   181
termSink.x=[211 189 330 176]/346*xL;
termSink.y=[173 184 192 218]/448*yL;
termSink.z=[ 75 126 133 181]/319*zL;
termSink.p0=133*10*[1 1 1 1];
termSink.perm=1e-6*[1 1 1 1];
    
% BigBrain Arterier:  dim(346 x 448 x 319)
% % %    185   180    98
% % %    171   184   144
termSrc.x=[185 171]/346*xL;
termSrc.y=[180 184]/448*yL;
termSrc.z=[ 98 144]/319*zL;
termSrc.p0=133*85*[1 1];
termSrc.perm=1e-6*[1 1];

%Viscosity
visc0 =0.004;
visc = visc0*ones(nx,ny,nz,2); % 4*1e-3*ones(nx,ny,nz,2); % the value is 3e-6 here..

% ### Voxel faces: (trans not perm) ...

%Permeability in x-direction
permX = 2e-10*ones(nx+1,ny,nz,2)/visc0;
%Permeability in y-direction
permY = 2e-10*ones(nx,ny+1,nz,2)/visc0;
%Permeability in z-direction
permZ = 2e-10*ones(nx,ny,nz+1,2)/visc0;

% % permX(:,:,:,2)=2e-10/visc0;
% % permY(:,:,:,2)=2e-10/visc0;
% % permZ(:,:,:,2)=2e-10/visc0;

%Inter-compartment conductivities
%This is a dimensionless quantity: Permeability [m^2] divided by
%interfacial area in a one one unit volume
%permQ = 2e-10*ones(nx,ny,nz);
permQ = 1e-8*ones(nx,ny,nz)/visc0;

% Terminal sources:
options.termSink.i=1+floor(nx*termSink.x/xL);
options.termSink.j=1+floor(ny*termSink.y/yL);
options.termSink.k=1+floor(nz*termSink.z/zL);
options.termSink.p0=termSink.p0;
options.termSink.perm=termSink.perm/visc0;
clear termSink;
options.termSrc.i=1+floor(nx*termSrc.x/xL);
options.termSrc.j=1+floor(ny*termSrc.y/yL);
options.termSrc.k=1+floor(nz*termSrc.z/zL);
options.termSrc.p0=termSrc.p0;
options.termSrc.perm=termSrc.perm/visc0;
clear termSrc;

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

% note that this overides any corresponding prs cnd
options.noFlowWestQ1 = 1;
options.noFlowEastQ1 = 1;
options.noFlowSouthQ1 = 1;
options.noFlowNorthQ1 = 1;
options.noFlowDownQ1 = 1; %close this to test terminal sources
options.noFlowUpQ1 = 1;
options.noFlowDownQ2 = 1; %close this to test terminal sinks
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

[state]=runFlowSolver3D2QSort(options);

