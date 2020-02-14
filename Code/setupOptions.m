%% create options for runFlowSolver3D2Q
clear options

%prm=rmfield(prm,'etafunc');
% here is the mask:

% time
options.time=0:0.1:60;
%options.time=0:0.1:60;
%options.time=0:0.001:60;
%Discretization
nx = 62; %prm.dim(1);
ny = 66;  %prm.dim(2);
nz = 20; %prm.dim(3); %3;

%Domain
xL = 4 * prm.h(1); % mm
yL = 4 * prm.h(2);
zL = 3 * prm.h(3);



%Viscosity
visc = prm.scaling*prm.tissue.visc*ones(nx,ny,nz,2); % 4*1e-3*ones(nx,ny,nz,2); % the value is 3e-6 here..

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
