
function drawbitbox3(numbits, boxcolors)

%set limits
xlim([0 3]);
ylim([0 numbits+2]);
yl = max(ylim);
xl = max(xlim);
%draw boxes
for b=0:numbits-10
    rectangle('position',[2,b,1,1],'edgecolor','k','facecolor',[1 1 1]*boxcolors(b+1));
end
for b=numbits-8:numbits-8+7
    rectangle('position',[2,b,1,1],'edgecolor','k','facecolor',[1 1 1]*boxcolors(b));
end
rectangle('position',[2,numbits+1,1,1],'edgecolor','k','facecolor',[1 1 1]*boxcolors(end));

text(-0.25,numbits+1+.4,'S','interpreter','latex','fontsize',7,'horizontalalignment','center');
text(-0.25,numbits-4+.4,'E','interpreter','latex','fontsize',7,'horizontalalignment','center');
text(-0.25,(23-(32-numbits))/2+.4,'F','interpreter','latex','fontsize',7,'horizontalalignment','center');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
axis off