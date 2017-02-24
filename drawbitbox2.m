
function drawbitbox2(numbits, bitblack)

%set limits
xlim([0 3]);
ylim([0 numbits+2]);
yl = max(ylim);
xl = max(xlim);

%texts dots
text(2.2,0+.5,'...','interpreter','tex','fontsize',8,'horizontalalignment','center','rotation',-90);

%draw boxes
txnum = numbits-9;
for b=1:numbits-10
    txnum = txnum-1;
    %rectangle('position',[2,b,1,1],'edgecolor','k','facecolor','w');
    text(2.2,b+.5,['F' num2str(txnum)],'interpreter','tex','fontsize',8,'horizontalalignment','center');
end


txnum = 9;
for b=numbits-8:numbits-8+7
    txnum = txnum-1;
    %rectangle('position',[2,b,1,1],'edgecolor','k','facecolor','w');
    text(3,b+.5,['E' num2str(txnum)],'interpreter',...
        'tex','fontsize',8,'horizontalalignment','right');
end
%rectangle('position',[2,numbits+1,1,1],'edgecolor','k','facecolor','w');
text(3,numbits+1.5,'Sign','interpreter','tex','fontsize',8,'horizontalalignment','right');

%text(-0.25,numbits+1+.4,'s','interpreter','latex','fontsize',10,'horizontalalignment','center');
%text(-0.25,numbits-4+.4,'e','interpreter','latex','fontsize',10,'horizontalalignment','center');
%text(-0.25,(23-(32-numbits))/2+.4,'f','interpreter','latex','fontsize',8,'horizontalalignment','center');



%fill bit black box
if bitblack>0
if bitblack==numbits
    bb=numbits+1;
elseif bitblack>numbits-9
    bb=bitblack;
else
    bb=bitblack-1;
end
rectangle('position',[2,bb,1,1],'edgecolor','none','facecolor','k');
end

%remove ticks
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);

axis off