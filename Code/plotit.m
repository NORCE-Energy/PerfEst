
swat = load('ADSO_TL50/swat');
poly = load('ADSO_TL50/poly');
eclStruct = importdata('ADSO_TL50/BASE120T005.RSM',' ',9);
ecl = eclStruct.data;

nFrames = size(ecl,1);
mov(1:nFrames) = struct('cdata', [],'colormap', []);
set(gca,'nextplot','replacechildren');

for index=[1:1:size(ecl,1)]
      i = cast(ecl(index,1),'int16');
      %%if (i > 2200), break, end
      dd1 = load(sprintf('polymer_output/concentration-%03d.dat', i));
      dd2 = load(sprintf('polymer_output/saturation-%03d.dat', i));
      subplot(2,1,1), plot(dd1(:), 'r'), hold on, plot(poly(:,index),'k'), axis([0,120,-0.1,5.1]), hold off
      title(['Polymer concentration at ',int2str(i),' days'],'fontsize',12,'fontweight','b')
      subplot(2,1,2), plot(dd2(1:2:end), 'r'), hold on, plot(swat(:,index),'k'), axis([0,120,-0.1,1.1]), hold off
      title(['Water saturation at ',int2str(i),' days'],'fontsize',12,'fontweight','b')
      bccn40(index)=dd1(40);
      bccn80(index)=dd1(80);
      satW40(index)=dd2(79);
      satW80(index)=dd2(159);
      %%subplot(5,1,4), plot(dd2(:,4), 'g'), hold on, axis([0,100,-0.1, 500.0]), hold off
      %%subplot(5,1,5), plot(dd2(:,5), 'b'), hold on, plot(dd2(:,6), 'r'), hold on, plot(dd2(:,7), 'g'), hold on, axis([0,100,-0.1,1.1]), hold off
  drawnow
  pause(0.1)
  mov(index) = getframe(gcf);
end

subplot(2,2,1), plot(ecl(:,1),bccn40(:), 'r'), hold on,  plot(ecl(:,1),ecl(:,6),'k'), axis([0,i,-0.1,5.5]), hold off
title('Polymer concentration at position 1/3','fontsize',12,'fontweight','b')
leg1 = legend('opm::poly','eclipse');
subplot(2,2,2), plot(ecl(:,1),bccn80(:), 'r'), hold on,  plot(ecl(:,1),ecl(:,7),'k'), axis([0,i,-0.1,5.5]), hold off
title('Polymer concentration at position 2/3','fontsize',12,'fontweight','b')
leg2 = legend('opm::poly','eclipse');
subplot(2,2,3), plot(ecl(:,1),satW40(:), 'r'), hold on,  plot(ecl(:,1),ecl(:,4),'k'), axis([0,i,-0.1,1.1]), hold off
title('Water saturation at position 1/3','fontsize',12,'fontweight','b')
leg3 = legend('opm::poly','eclipse');
subplot(2,2,4), plot(ecl(:,1),satW80(:), 'r'), hold on,  plot(ecl(:,1),ecl(:,5),'k'), axis([0,i,-0.1,1.1]), hold off
title('Water saturation at position 2/3','fontsize',12,'fontweight','b')
leg4 = legend('opm::poly','eclipse');

movie2avi(mov, 'noAds_noDead.avi', 'compression', 'None');
