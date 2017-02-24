
function bitinfoplot(p1,p2,p3,q2,q3,scy,bins)
    
p1(p1==0) = .000000001;
p2(p2==0) = .000000001;
p3(p3==0) = .000000001;


entr1 = sum(-p1.*log2(p1));
entr2 = sum(-p2.*log2(p2));
entr3 = sum(-p3.*log2(p3));

info = entr1-q2*entr2-q3*entr3;
if q2==0 || q3==0
    info=0;
end

hold on
plot(bins,p2,'color',[1 1 1]*0,'linestyle','-'); 
plot(bins,p3,'color',[1 1 1]*.5,'linestyle','-'); 
plot(bins,p1,'color',[1 1 1]*0,'linestyle',':'); 

ylim([0 scy*max(p1)]);

hold off

set(gca,'xtick',[]);
set(gca,'ytick',[]);

text(double(min(xlim)+.67*(max(xlim)-min(xlim))),double(.89*max(ylim)),['$I$=' num2str(info,'%1.3f')],'fontsize',10,'interpreter','latex');
