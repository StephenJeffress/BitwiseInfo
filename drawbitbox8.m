
function drawbitbox8(numbits, blckbx)

boxcolors = ones(1,numbits);
boxcolors(blckbx) = 0;

xlim([0 4]);
ylim([0 numbits+2]);
yl = max(ylim);
xl = max(xlim);
%draw boxes
for b=0:10
    rectangle('position',[2,1.7*b,2,1.7],'edgecolor','k','facecolor',[1 1 1]*boxcolors(b+1));
end
for b=11:numbits-10
    rectangle('position',[0,1.7*b-1.7*11,2,1.7],'edgecolor','k','facecolor',[1 1 1]*boxcolors(b+1));
end

for b=numbits-8:numbits-8+3
    rectangle('position',[2,1.7*b-1.7*11+.5,2,1.7],'edgecolor','k','facecolor',[1 1 1]*boxcolors(b));
end
for b=numbits-8+4:numbits-8+7
    rectangle('position',[0,1.7*b-1.7*15+.5,2,1.7],'edgecolor','k','facecolor',[1 1 1]*boxcolors(b));
end
rectangle('position',[2,numbits+.1,2,1.7],'edgecolor','k','facecolor',[1 1 1]*boxcolors(end));
text(-5.5,numbits+1+.4,'Sign','interpreter','tex','fontsize',8,'horizontalalignment','left');
%text(-2.5,numbits+1+.4,'Sign','interpreter','tex','fontsize',6,'horizontalalignment','center');
%text(-6.5,numbits-6+.4,'Exp.','interpreter','tex','fontsize',8,'horizontalalignment','left');
text(-3.5,numbits-9+.4,'Exp','interpreter','tex','fontsize',8,'horizontalalignment','left','rotation',90);
text(-3.5,numbits-27+.4,'Fraction','interpreter','tex','fontsize',8,'horizontalalignment','left','rotation',90);
% text(-0.25,(23-(32-numbits))/2+.4,'F','interpreter','latex','fontsize',7,'horizontalalignment','center');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
axis off

