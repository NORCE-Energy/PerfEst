function [func,grad] = funcGrad(logPermXYZ, pars)

options = pars.options;
nx = options.nx;
ny = options.ny;
nz = options.nz;
nn=nx*ny*nz;
options.permX(:,:,:,1)=reshape(exp(logPermXYZ(1:nn)),nx,ny,nz);
options.permY(:,:,:,1)=reshape(exp(logPermXYZ(nn+1:2*nn)),nx,ny,nz);
options.permZ(:,:,:,1)=reshape(exp(logPermXYZ(2*nn+1:3*nn)),nx,ny,nz);

options.porosityArt = reshape(logPermXYZ(3*nn+1:4*nn),nx,ny,nz);
disp('### funcGrad:  include porosity ###');

% if (max(logPermXYZ(3*nn+1:4*nn))>1.0 || min(logPermXYZ(3*nn+1:4*nn))<0.0000001)
%   func=1e10;
%   grad=zeros(4*nn,1);
%   return;
% end

[state]=runFlowSolver3D2QSort(options);

[gradPerm, gradPor, sampleMask, lamTau, lamSrc, lamPrs] = computGradAdj(options, state, pars.tofArtObserved, pars.alphaPermReg, 0);

d_tofArt=(state.tofArt-pars.tofArtObserved).*reshape(sampleMask,nx,ny,nz);
%func = 0.5*reshape(d_tofArt,1,[])*reshape(d_tofArt,1,[])';
porA=options.porosityArt(:);
gamma=100.0;
func=porA.*(porA<0)+(porA-1).*(porA>1);
func=0.5*gamma*func'*func + 0.5*reshape(d_tofArt,1,[])*reshape(d_tofArt,1,[])';

if (pars.alphaPermReg ~= 0)    %regularization terms
  func = func + (pars.alphaPermReg*norm(reshape(options.permX(:,:,:,1)-2e-10,1,[]).*sampleMask))^2 ...
              + (pars.alphaPermReg*norm(reshape(options.permY(:,:,:,1)-2e-10,1,[]).*sampleMask))^2 ...
              + (pars.alphaPermReg*norm(reshape(options.permZ(:,:,:,1)-2e-10,1,[]).*sampleMask))^2;
end

% % % if (0==1)
% % %   for j=1:10
% % %     climsTofErr = [-0.2 0.2];
% % %     subplot(10,10,10*(j-1)+7);
% % %     imagesc(d_tofArt(:,:,j),climsTofErr);
% % %     if j==1
% % %       c=colorbar;
% % %       c.Location='north';
% % %     end
% % %     ax = gca;
% % %     ax.XAxis.Visible = 'off';
% % %     ax.YAxis.Visible = 'off';
% % % 
% % %     climsLamTau = [-1.0e9 1.0e9];
% % %     subplot(10,10,10*(j-1)+8);
% % %     imagesc(lamTau(:,:,j),climsLamTau);
% % %     if j==1
% % %       c=colorbar;
% % %       c.Location='northoutside';
% % %     end
% % %     ax = gca;
% % %     ax.XAxis.Visible = 'off';
% % %     ax.YAxis.Visible = 'off';   
% % % 
% % %     climsLamSrc = [-30 30];
% % %     subplot(10,10,10*(j-1)+9);
% % %     imagesc(lamSrc(:,:,j),climsLamSrc);
% % %     if j==1
% % %       c=colorbar;
% % %       c.Location='northoutside';
% % %     end
% % %     ax = gca;
% % %     ax.XAxis.Visible = 'off';
% % %     ax.YAxis.Visible = 'off';   
% % % 
% % %     climsLamPrs = [-10e11 -5e11];
% % %     subplot(10,10,10*(j-1)+10);
% % %     imagesc(lamPrs(:,:,j,1),climsLamPrs);
% % %     if j==1
% % %       c=colorbar;
% % %       c.Location='northoutside';
% % %     end
% % %     ax = gca;
% % %     ax.XAxis.Visible = 'off';
% % %     ax.YAxis.Visible = 'off';   
% % % 
% % % 
% % %     set(gca,'YDir','normal')
% % %   end
% % % 
% % %   drawnow,pause(1)
% % % end

% perm
gradPermX = reshape(reshape(gradPerm(:,:,:,1,1),1,[]),nx,ny,nz);
gradPermY = reshape(reshape(gradPerm(:,:,:,1,2),1,[]),nx,ny,nz);
gradPermZ = reshape(reshape(gradPerm(:,:,:,1,3),1,[]),nx,ny,nz);
  
% log perm ...
gradLogPermX = gradPermX .* options.permX(:,:,:,1);
gradLogPermY = gradPermY .* options.permY(:,:,:,1);
gradLogPermZ = gradPermZ .* options.permZ(:,:,:,1);

%limit gradents
gMax=10;
gradLogPermX = min(gradLogPermX, gMax);
gradLogPermX = max(gradLogPermX, -gMax);
gradLogPermY = min(gradLogPermY, gMax);
gradLogPermY = max(gradLogPermY, -gMax);
  
gradLogPermXAbsSq = reshape(gradLogPermX,1,[])*reshape(gradLogPermX,1,[])';
gradLogPermYAbsSq = reshape(gradLogPermY,1,[])*reshape(gradLogPermY,1,[])';
gradLogPermZAbsSq = reshape(gradLogPermZ,1,[])*reshape(gradLogPermZ,1,[])';
gradLogPermAbsSq=gradLogPermXAbsSq+gradLogPermYAbsSq+gradLogPermZAbsSq+1.0e-25;
  
grad(1:nn) = reshape(gradLogPermX,1,[]);
grad(nn+1:2*nn) = reshape(gradLogPermY,1,[]);
grad(2*nn+1:3*nn) = reshape(gradLogPermZ,1,[]);

%gradPor=0.0001*gradPor;

%%% Denne gir mer "struktur", men ikke vesentlig bedre tilpassning ...
% alphaPor=0.05;
% gradPor = min(gradPor, 0.99*options.porosityArt);
% gradPor = alphaPor*max(gradPor, -0.99*(1-options.porosityArt));

%grad(1:3*nn)=0;
%gradPor=gradPor + gamma*(porA.*(porA<0)+(porA-1).*(porA>1));

%gradPor = zeros(nn,1);
%%% Ok tilpassning gitt "rett" H0_phi, men mindre tydlig struktur ...
grad(3*nn+1:4*nn) = reshape(gradPor,1,[]) + (gamma*(porA.*(porA<0)+(porA-1).*(porA>1)))';
disp('### funcGrad:  include pososity ###');

grad=grad';
  
msg=['          f(x): ',num2str(func), ...
     ' gradPor: ', num2str(norm(reshape(gradPor,1,[]))), ...
     ' gradLogPermX: ', num2str(sqrt(gradLogPermXAbsSq)), ...
     ' gradLogPermY: ', num2str(sqrt(gradLogPermYAbsSq)), ...
     ' gradLogPermAbs: ',num2str(sqrt(gradLogPermAbsSq))];
disp(msg);

plotfig = false;
if plotfig
    figure('Name','frog','NumberTitle','off');
    subplot(2,5,1);
    imagesc(d_tofArt);
    title('dTofArt')
    c=colorbar;
    subplot(2,5,2);
    imagesc(lamTau);
    title('lamTauArt')
    c=colorbar;
    subplot(2,5,3);
    imagesc(lamSrc(:,:,1,1));
    title('lamSrcArt')
    c=colorbar;
    subplot(2,5,4);
    imagesc(lamPrs(:,:,1,1));
    title('lamPrsArt')
    c=colorbar;
    subplot(2,5,5);
    imagesc(gradPor(:,:,1));
    title('gradPorArt')
    c=colorbar;
    subplot(2,5,6);
    imagesc(options.porosityArt(:,:,1));
    title('porosityArt')
    c=colorbar;
    subplot(2,5,7);
    imagesc(gradLogPermX);
    title('gradLogPermXArt')
    c=colorbar;
    subplot(2,5,8);
    imagesc(gradLogPermY);
    title('gradLogPermYArt')
    c=colorbar;
    subplot(2,5,9);
    imagesc(log(options.permX(:,:,1,1)));
    title('permXArtLog')
    c=colorbar;
    subplot(2,5,10);
    imagesc(log(options.permY(:,:,1,1)));
    title('permYArtLog')
    c=colorbar;
    
    lamSrcArt=lamSrc(:,:,1,1);
    lamPrsArt=lamPrs(:,:,1,1);
    permX=options.permX(:,:,1,1);
    permY=options.permY(:,:,1,1);
    porosityArt=options.porosityArt(:,:,1);
    
    % % for j=1:nz
    % %   climsPerm = [-11 -5];
    % %   subplot(nz,10,10*(j-1)+1);
    % %   imagesc(log10(options.permX(:,:,j,1)),climsPerm);
    % %   if j==1
    % %     c=colorbar;
    % %     c.Location='northoutside';
    % %   end
    % %   ax = gca;
    % %   ax.XAxis.Visible = 'off';
    % %   ax.YAxis.Visible = 'off';
    % %   subplot(nz,10,10*(j-1)+2);
    % %   imagesc(log10(options.permY(:,:,j,1)),climsPerm);
    % %   ax = gca;
    % %   ax.XAxis.Visible = 'off';
    % %   ax.YAxis.Visible = 'off';
    % %   subplot(nz,10,10*(j-1)+3);
    % %   imagesc(log10(options.permZ(:,:,j,1)),climsPerm);
    % %   ax = gca;
    % %   ax.XAxis.Visible = 'off';
    % %   ax.YAxis.Visible = 'off';
    % %
    % %   climsGrad = [-0.5 0.5];
    % %   subplot(nz,10,10*(j-1)+4);
    % %   imagesc(gradLogPermX(:,:,j),climsGrad);
    % %   if j==1
    % %     c=colorbar;
    % %     c.Location='northoutside';
    % %   end
    % %   ax = gca;
    % %   ax.XAxis.Visible = 'off';
    % %   ax.YAxis.Visible = 'off';
    % %   subplot(nz,10,10*(j-1)+5);
    % %   imagesc(gradLogPermY(:,:,j),climsGrad);
    % %   ax = gca;
    % %   ax.XAxis.Visible = 'off';
    % %   ax.YAxis.Visible = 'off';
    % %   subplot(nz,10,10*(j-1)+6);
    % %   imagesc(gradLogPermZ(:,:,j),climsGrad);
    % %   ax = gca;
    % %   ax.XAxis.Visible = 'off';
    % %   ax.YAxis.Visible = 'off';
    % %
    % %   climsTofErr = [-0.2 0.2];
    % %   subplot(nz,10,10*(j-1)+7);
    % %   imagesc(d_tofArt(:,:,j),climsTofErr);
    % %   if j==1
    % %     c=colorbar;
    % %     c.Location='north';
    % %   end
    % %   ax = gca;
    % %   ax.XAxis.Visible = 'off';
    % %   ax.YAxis.Visible = 'off';
    % %
    % %   climsLamTau = [-1.0e9 1.0e9];
    % %   subplot(nz,10,10*(j-1)+8);
    % %   imagesc(lamTau(:,:,j),climsLamTau);
    % %   if j==1
    % %     c=colorbar;
    % %     c.Location='northoutside';
    % %   end
    % %   ax = gca;
    % %   ax.XAxis.Visible = 'off';
    % %   ax.YAxis.Visible = 'off';
    % %
    % %   climsLamSrc = [-30 30];
    % %   subplot(nz,10,10*(j-1)+9);
    % %   imagesc(lamSrc(:,:,j),climsLamSrc);
    % %   if j==1
    % %     c=colorbar;
    % %     c.Location='northoutside';
    % %   end
    % %   ax = gca;
    % %   ax.XAxis.Visible = 'off';
    % %   ax.YAxis.Visible = 'off';
    % %
    % %   climsLamPrs = [1e13 3e13];
    % %   subplot(nz,10,10*(j-1)+10);
    % %   imagesc(lamPrs(:,:,j,1),climsLamPrs);
    % %   if j==1
    % %     c=colorbar;
    % %     c.Location='northoutside';
    % %   end
    % %   ax = gca;
    % %   ax.XAxis.Visible = 'off';
    % %   ax.YAxis.Visible = 'off';
    % %
    % %   set(gca,'YDir','normal')
    % %   %colormap(hot)
    % % end
    
    drawnow,pause(1)
    
    closereq;
end

% figure('Name','perm X Y Z','NumberTitle','off')
% semilogy(reshape(options.permX(5,5,:,1),1,[]),'DisplayName','permX');
% hold on
% semilogy(reshape(options.permY(5,5,:,1),1,[]),'DisplayName','permY');
% semilogy(reshape(options.permZ(5,5,:,1),1,[]),'DisplayName','permZ');
% legend('show');
% hold off

% permX = options.permX(:,:,1:10,1)
% permY = options.permY(:,:,1:10,1)
% permZ = options.permZ(:,:,1:10,1)
      
end


  
  
  
  
  
