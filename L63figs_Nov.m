

%%

%initialize vars
clc; clear; clf; %close all;

tic

savefigs = 1;
dockfigs = 0;
%savedir = '~/public_html/';
%savedir = '/home/jeffress/Documents/Publications/BitwiseInfo/19Nov2016/';
savedir = '/home/jeffress/Documents/Publications/BitwiseInfo/10Jan2017/';

vartype = 'single';
%vartype = 'double';


K= 3;
a = ([10, 28, 8/3]);  %sigma, rho, beta
%a = [20, 600, 3];  %sigma, rho, beta


dt =             0.01;
steps =          1e8;    %50000000;
psize =          0.3;
TNres =          5     ;
bwid = psize;
pevery =            0.001;


if strcmp(vartype,'single');
    dt = single(    dt    );
    %steps = single(  steps  );
    psize = single(  psize     );
    TNres = single(  TNres     );
    a = single(a);
    K = single(K);
end



Plotlength = round(25/dt);
bits2=15;
bits3=12;
bits4=9;

Numbits = 32-9;

%allocate memory
dX = zeros(K,1,vartype);
K1=zeros(K,1,vartype);
K2=zeros(K,1,vartype);
K3=zeros(K,1,vartype);
K4=zeros(K,1,vartype);
X=zeros(K,steps,vartype);
Xp=zeros(K,Plotlength,vartype);
T = zeros(1,steps,vartype);

% %randomly initialize
% for i=1:K
%     X(i,1)=rand;
% end
% %run for 2000 to get a better initial condition
% for i=1:5000;
%    X(:,i+1)=X(:,i)+RungeLorenz63(X(:,i),dt,dX,K1,K2,K3,K4,a);
% end

load x0.mat;
X(:,1) = x;
%X(:,1) = X(:,4000);
Xp(:,1)=X(:,1)+psize*ones(3,1);


%run full experiment
nup=single(0);
for i=1:steps-1;
   Dx = RungeLorenz63(X(:,i),dt,dX,K1,K2,K3,K4,a);
   %Dx = round(Dx)
   X(:,i+1)=X(:,i)+Dx;
   
   if i<Plotlength
       Dx = RungeLorenz63(Xp(:,i),dt,dX,K1,K2,K3,K4,a);
       %Dx = round(Dx)
       Xp(:,i+1)=Xp(:,i)+Dx;
   end
   
   if i<steps-3
   if X(1,i+1)==0 || X(2,i+1)==0 || X(3,i+1)==0 
       disp('why');
       X(:,i+1)=X(:,i-1)+randn(3,1);
   end
   end
   
%    X(1,i+1) = trimbits2(X(1,i+1),3);
%    X(2,i+1) = trimbits2(X(2,i+1),3);
%    X(3,i+1) = trimbits2(X(3,i+1),3);
%    
   T(i+1)=T(i)+dt;
   
   if i/steps > nup
       nup = nup+pevery
   end
   
end
%Xp=Xp./tstep;

% 
% clf
% plot3(X(1,:),X(2,:),X(3,:))

%
%

%%


Xb = getbits(X(1,:));
Yb = getbits(X(2,:));
Zb = getbits(X(3,:));

%%

%save l63longrun.mat X Xb Yb Zb -v7.3


%% FIG 1 data

%get all vars needed for plots


X1colonPln = X(1,1:Plotlength);
X2colonPln = X(2,1:Plotlength);
X3colonPln = X(3,1:Plotlength);
%X1pcolonPln = Xp(1,1:Plotlength);
tmpn=round(3.5*Plotlength);
X1colonPln2 = X(1,1:tmpn);
X2colonPln2 = X(2,1:tmpn);
X3colonPln2 = X(3,1:tmpn);

T1colonPln=T(1:Plotlength);
blcks = 100;
bln = round(linspace(1,Plotlength,blcks));
Bitimage = Xb(bln,:);


bwid=.3;

tmsend = 5.25e7;

tv = X(1,1:tmsend);
Bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
%Bins = [-fliplr(bwid:bwid:abs(min(tv))+bwid) (0:bwid50000000:abs(max(tv))+bwid)];
pdf = single(hist(tv,Bins));
X1pdf = pdf./sum(pdf);


bnum=32;
%fdn = 1;
tv = X(1,1:tmsend);
cv = squeeze(Xb(1:length(tv),bnum));

Psign1 = sum(cv==true)/length(cv);
Xsign1pdf = hist(tv(cv==true),Bins);
Xsign1pdf = Xsign1pdf/sum(Xsign1pdf);

Psign0 = sum(cv==false)/length(cv);
Xsign0pdf = hist(tv(cv==false),Bins);
Xsign0pdf = Xsign0pdf/sum(Xsign0pdf);
Xsign1pdf(Xsign1pdf==0) = .000000001;
Xsign0pdf(Xsign0pdf==0) = .000000001;


%% Fig 2 data


tmsend = 5.25e7;

t1=1;t2=round(.2/dt);t3=round(1/dt);
tv = X(1,1:tmsend);


bnum=32;
cv = squeeze(Xb(1:length(tv),bnum));
p1baseline = single(hist(tv,Bins));
p1baseline = p1baseline./sum(p1baseline);
q2b32_t1 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b32_t1 = p2/sum(p2);
q3b32_t1 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b32_t1 = p3/sum(p3);


tv = X(1,t2:tmsend);
cv = squeeze(Xb(1:length(tv),bnum));
%p1 = single(hist(tv,Bins));
%p1 = p1./sum(p1);
q2b32_t2 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b32_t2 = p2/sum(p2);
q3b32_t2 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b32_t2 = p3/sum(p3);


tv = X(1,t3:tmsend);
cv = squeeze(Xb(1:length(tv),bnum));
%p1 = single(hist(tv,Bins));
%p1 = p1./sum(p1);
q2b32_t3 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b32_t3 = p2/sum(p2);
q3b32_t3 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b32_t3 = p3/sum(p3);



bnum=31;
tv = X(1,t1:tmsend);
cv = squeeze(Xb(1:length(tv),bnum));
%p1 = single(hist(tv,Bins));
%p1 = p1./sum(p1);
q2b31_t1 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b31_t1 = p2/sum(p2);
q3b31_t1 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b31_t1 = p3/sum(p3);

tv = X(1,t2:tmsend);
cv = squeeze(Xb(1:length(tv),bnum));
q2b31_t2 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b31_t2 = p2/sum(p2);
q3b31_t2 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b31_t2 = p3/sum(p3);

tv = X(1,t3:tmsend);
cv = squeeze(Xb(1:length(tv),bnum));
q2b31_t3 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b31_t3 = p2/sum(p2);
q3b31_t3 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b31_t3 = p3/sum(p3);




bnum=24;
tv = X(1,t1:end);
cv = squeeze(Xb(1:length(tv),bnum));
%p1 = single(hist(tv,Bins));
%p1 = p1./sum(p1);
q2b24_t1 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b24_t1 = p2/sum(p2);
q3b24_t1 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b24_t1 = p3/sum(p3);

tv = X(1,t2:end);
cv = squeeze(Xb(1:length(tv),bnum));
q2b24_t2 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b24_t2 = p2/sum(p2);
q3b24_t2 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b24_t2 = p3/sum(p3);

tv = X(1,t3:end);
cv = squeeze(Xb(1:length(tv),bnum));
q2b24_t3 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b24_t3 = p2/sum(p2);
q3b24_t3 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b24_t3 = p3/sum(p3);


bnum=23;
tv = X(1,t1:end);
cv = squeeze(Xb(1:length(tv),bnum));
%p1 = single(hist(tv,Bins));
%p1 = p1./sum(p1);
q2b23_t1 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b23_t1 = p2/sum(p2);
q3b23_t1 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b23_t1 = p3/sum(p3);

tv = X(1,t2:end);
cv = squeeze(Xb(1:length(tv),bnum));
q2b23_t2 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b23_t2 = p2/sum(p2);
q3b23_t2 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b23_t2 = p3/sum(p3);

tv = X(1,t3:end);
cv = squeeze(Xb(1:length(tv),bnum));
q2b23_t3 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b23_t3 = p2/sum(p2);
q3b23_t3 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b23_t3 = p3/sum(p3);


bnum=20;
tv = X(1,t1:end);
cv = squeeze(Xb(1:length(tv),bnum));
%p1 = single(hist(tv,Bins));
%p1 = p1./sum(p1);
q2b20_t1 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b20_t1 = p2/sum(p2);
q3b20_t1 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b20_t1 = p3/sum(p3);

tv = X(1,t2+2:end);
cv = squeeze(Xb(1:length(tv),bnum));
q2b20_t2 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b20_t2 = p2/sum(p2);
q3b20_t2 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b20_t2 = p3/sum(p3);

tv = X(1,t3:end);
cv = squeeze(Xb(1:length(tv),bnum));
q2b20_t3 = sum(cv==true)/length(cv);
p2 = hist(tv(cv==true),Bins);
p2b20_t3 = p2/sum(p2);
q3b20_t3 = sum(cv==false)/length(cv);
p3 = hist(tv(cv==false),Bins);
p3b20_t3 = p3/sum(p3);



%% Fig 3 data

Numbits=22;

timeTN = 10.^linspace(log10(dt),log10(200),45);
TNs = ceil(timeTN./dt);
TNs(1:11) = 1:11;

%Tinfo = zeros(Numbits,length(TNs));
Tinfox = zeros(Numbits,length(TNs));

%for locK = 1:3

for Tn = 1:length(TNs)
    tn = TNs(Tn); 
for bnum=1:Numbits
    [bnum Tn/length(TNs)]
    tv = X(1,tn:end);
    Bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
    cv = squeeze(Xb(1:length(tv),32-Numbits+bnum));
    tmpinfo = getinfo(tv,cv,Bins);
    %Tinfo(bnum,Tn) = Tinfo(bnum,Tn)+tmpinfo;
    Tinfox(bnum,Tn) = Tinfox(bnum,Tn)+tmpinfo;
end
end

% for Tn = 1:length(TNs)
%     [locK 2/3 Tn/length(TNs)]
%     tn = TNs(Tn); 
% for bnum=1:Numbits
%     tv = X(locK,tn:end);
%     Bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
%     cv = squeeze(Yb(1:length(tv),32-Numbits+bnum));
%     Tinfo(bnum,Tn) = Tinfo(bnum,Tn)+getinfo(tv,cv,Bins);
% end
% end
% for Tn = 1:length(TNs)
%     [locK 3/3 Tn/length(TNs)]
%     tn = TNs(Tn); 
% for bnum=1:Numbits
%     tv = X(locK,tn:end);
%     Bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
%     cv = squeeze(Zb(1:length(tv),32-Numbits+bnum));
%     Tinfo(bnum,Tn) = Tinfo(bnum,Tn)+getinfo(tv,cv,Bins);
% end
% end
% end


%%


%save vars
save l63plotvars.mat X1colonPln X2colonPln X3colonPln ...
    X1colonPln2 X2colonPln2 X3colonPln2 X1pcolonPln ...
    T1colonPln Bitimage Bins X1pdf Psign1 Xsign1pdf ...
    Psign0 Xsign0pdf Plotlength dt steps psize TNres bwid ...
p1baseline  p2b32_t1 p3b32_t1 q2b32_t1 q3b32_t1 ...
    p2b32_t2 p3b32_t2 q2b32_t2 q3b32_t2 ...
    p2b32_t3 p3b32_t3 q2b32_t3 q3b32_t3 ...
    p2b31_t1 p3b31_t1 q2b31_t1 q3b31_t1 ...
    p2b31_t2 p3b31_t2 q2b31_t2 q3b31_t2 ...
    p2b31_t3 p3b31_t3 q2b31_t3 q3b31_t3 ...
    p2b24_t1 p3b24_t1 q2b24_t1 q3b24_t1 ...
    p2b24_t2 p3b24_t2 q2b24_t2 q3b24_t2 ...
    p2b24_t3 p3b24_t3 q2b24_t3 q3b24_t3 ...
    p2b23_t1 p3b23_t1 q2b23_t1 q3b23_t1 ...
    p2b23_t2 p3b23_t2 q2b23_t2 q3b23_t2 ...
    p2b23_t3 p3b23_t3 q2b23_t3 q3b23_t3 ...
    p2b20_t1 p3b20_t1 q2b20_t1 q3b20_t1 ...
    p2b20_t2 p3b20_t2 q2b20_t2 q3b20_t2 ...
    p2b20_t3 p3b20_t3 q2b20_t3 q3b20_t3 ...
    Tinfo Tinfox timeTN TNs Diffpf Diff2f ... 
    Diff3f Diff4f pdfXp pdfX2p pdfX3p pdfX4p Bins2...
    X Xp X2p X3p X4p fnum bits2 bits3 bits4 attX4p ...
    attX3p attX2p attX attXp;

%%


clear all
load l63plotvars.mat
%load TINFO.mat

dockfigs=0;
savefigs=1;
Numbits = 32;
%savedir = '/home/jeffress/Documents/Publications/BitwiseInfo/19Nov2016/';
savedir = '/home/jeffress/Documents/Publications/BitwiseInfo/10Jan2017/';


%% Fig 4

close all

set(0,'defaultaxesfontsize',10);
set(0,'defaulttextfontsize',10);

figure(1)
if dockfigs
    set(gcf,'windowstyle','docked');
else
    set(gcf,'visible','off');
end
clf

Plotlength = round(15/dt);

%Plotlength = length(X2p(1,:));

bxlf = .120; %.075
bxbt = .12;
bxwd = .20;
bxht = .26;
bxpw = .012;
bxpu = .02;


%Timeseries
axes('position',[bxlf bxbt+2*bxht+2*bxpu bxwd bxht])
hold on
plot(dt*(1:Plotlength),X(1,1:Plotlength),'k'); 
plot(dt*(1:Plotlength),X2p(1,1:Plotlength),'color',[1 1 1]*.5,'linestyle','-'); hold off;
hold off
ylim([-1 1]*45)
xlim([0 1]*15)
box on
%grid on
set(gca,'xtick',[0:5:15]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[0]);
set(gca,'yticklabel',[]);
set(gca,'xaxislocation','top');
xlabel('Timeseries','fontsize',10)
line(0+[.7 2],[1 1]*30.5,'color','k');
text(2.25,31,'$x(t)$','interpreter','latex','horizontalalignment','left','fontsize',10);
line(5+[.7 2],[1 1]*30.5,'color',[1 1 1]*.5);
text(7.25,31,'$x''_{12\:bit}(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');



axes('position',[bxlf bxbt+bxht+bxpu bxwd bxht])
hold on
plot(dt*(1:Plotlength),X(1,1:Plotlength),'k'); 
plot(dt*(1:Plotlength),X3p(1,1:Plotlength),'color',[1 1 1]*.5,'linestyle','-'); hold off;
hold off
ylim([-1 1]*45)
xlim([0 1]*15)
box on
%grid on
set(gca,'xtick',[0:5:10]);
set(gca,'ytick',[0]);
%set(gca,'ytick',[-20 -10 0 10 20]);
set(gca,'yticklabel',[]);
%set(gca,'xaxislocation','top');
set(gca,'xticklabel',[]);%grid on
%line(0+[.7 2],[1 1]*24.5,'color','k');
%text(2.25,25,'$x(t)$','interpreter','latex','horizontalalignment','left','fontsize',10);
line(0+[.7 2],[1 1]*30.5,'color',[1 1 1]*.5);
text(2.25,31,'$x''_{14\:bit}(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');



axes('position',[bxlf bxbt bxwd bxht])
box on
set(gca,'xtick',[]);set(gca,'ytick',[]);
hold on
plot(dt*(1:Plotlength),X(1,1:Plotlength),'k'); 
plot(dt*(1:Plotlength),X4p(1,1:Plotlength),'color',[1 1 1]*.5,'linestyle','-'); hold off;
hold off
ylim([-1 1]*45)
xlim([0 1]*15)
box on
%grid on
set(gca,'xtick',[0:5:15]);
set(gca,'ytick',[0]);
set(gca,'yticklabel',[]);
set(gca,'xaxislocation','bottom');
%set(gca,'xticklabel',[]);%grid on
h = xlabel('Time (mtu)','fontsize',10);
%xlbbot = get(h,'position');
%line(0+[.7 2],[1 1]*24.5,'color','k');
%text(2.25,25,'$x(t)$','interpreter','latex','horizontalalignment','left','fontsize',10);
line(0+[.7 2],[1 1]*30.5,'color',[1 1 1]*.5);
text(2.25,31,'$x''_{16\:bit}(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');






%Atrractors
atst = 12000;
atfn = 20000;
axes('position',[bxlf+1*bxwd+1*bxpw bxbt+2*bxht+2*bxpu bxwd bxht])
plot(attX2p(1,atst:atfn),attX2p(3,atst:atfn),'color',[1 1 1]*.5,'linewidth',.5);
xlim([-1 1]*22)
ylim([0 1]*55)
set(gca,'ytick',[]);
set(gca,'xtick',[-15 0 15]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
set(gca,'xaxislocation','top');
xlabel('Attractor','fontsize',10,'interpreter','tex');
line(min(xlim)+10+[1.5 4.5],[1 1]*.835*max(ylim),'color',[1 1 1]*.5);
text(min(xlim)+15.25,.835*max(ylim),'${\bf x}''_{12\:bit}(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');
box on

axes('position',[bxlf+1*bxwd+1*bxpw bxbt+bxht+bxpu bxwd bxht])
plot(attX3p(1,atst:atfn),attX3p(3,atst:atfn),'color',[1 1 1]*.5,'linewidth',.5);
xlim([-1 1]*22)
ylim([0 1]*55)
set(gca,'xtick',[-15 0 15]);
set(gca,'yticklabel',[]); 
set(gca,'xticklabel',[]);
line(min(xlim)+10+[1.5 4.5],[1 1]*.835*max(ylim),'color',[1 1 1]*.5);
text(min(xlim)+15.25,.835*max(ylim),'${\bf x}''_{14\:bit}(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');
box on

axes('position',[bxlf+1*bxwd+1*bxpw bxbt bxwd bxht])
plot(attX4p(1,atst:atfn),attX4p(3,atst:atfn),'color',[1 1 1]*.5,'linewidth',.5);
xlim([-1 1]*22)
ylim([0 1]*55)
set(gca,'xtick',[-15 0 15]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',{'-15' '0' '15'});
xlabel('$x$','fontsize',10,'interpreter','latex');
line(min(xlim)+10+[1.5 4.5],[1 1]*.835*max(ylim),'color',[1 1 1]*.5);
text(min(xlim)+15.25,.835*max(ylim),'${\bf x}''_{16\:bit}(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');
box on







%PDFS
axes('position',[bxlf+2*bxwd+2*bxpw bxbt+2*bxht+2*bxpu bxwd bxht])
hold on
plot(Bins2,pdfXp,'k');
plot(Bins2,pdfX2p,'color',[1 1 1]*.5,'linestyle','-');
xlim([-1 1]*22)
ylim([0 .02]);
set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',{'-15' '0' '15'}); 
set(gca,'ytick',[.2 .5 .8]*max(ylim));set(gca,'yticklabel',[]); 
hold off
box on
set(gca,'xaxislocation','top');
set(gca,'xticklabel',[]);
xlabel('Climatology','fontsize',10,'interpreter','tex');
ylim([0 .024])
line(min(xlim)+0+[1.5 4.5],[1 1]*.835*max(ylim),'color',[1 1 1]*0);
text(min(xlim)+5.25,.835*max(ylim),'$p(x)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');
line(min(xlim)+13+[1.5 4.5],[1 1]*.835*max(ylim),'color',[1 1 1]*.5);
text(min(xlim)+18.25,.835*max(ylim),'$p(x''_{12\:bit})$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');



axes('position',[bxlf+2*bxwd+2*bxpw bxbt+bxht+bxpu bxwd bxht])
hold on
plot(Bins2,pdfXp,'k');
plot(Bins2,pdfX3p,'color',[1 1 1]*.5,'linestyle','-');
xlim([-1 1]*22)
ylim([0 .02]);
set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',{'-15' '0' '15'}); 
set(gca,'ytick',[.2 .5 .8]*max(ylim));set(gca,'yticklabel',[]);
ylim([0 .024])
hold off
box on
line(min(xlim)+0+[1.5 4.5],[1 1]*.835*max(ylim),'color',[1 1 1]*.5);
text(min(xlim)+5.25,.835*max(ylim),'$p(x''_{14\:bit})$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');
axes('position',[bxlf+2*bxwd+2*bxpw bxbt bxwd bxht])
hold on
plot(Bins2,pdfXp,'k');
plot(Bins2,pdfX4p,'color',[1 1 1]*.5,'linestyle','-');
xlim([-1 1]*22)
ylim([0 .024])
set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',{'-15' '0' '15'}); 
set(gca,'ytick',[.2 .5 .8]*max(ylim));set(gca,'yticklabel',[]);
hold off
box on
xlabel('$x$','fontsize',10,'interpreter','latex');
line(min(xlim)+0+[1.5 4.5],[1 1]*.835*max(ylim),'color',[1 1 1]*.5);
text(min(xlim)+5.25,.835*max(ylim),'$p(x''_{16\:bit})$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');




%ERROR GROWTH
fnum=3e5;
axes('position',[bxlf+3*bxwd+3*bxpw bxbt+2*bxht+2*bxpu bxwd bxht])
semilogx((1:length(Diffpf))*dt,Diffpf/(fnum),'color','k','linestyle','-','linewidth',1); hold on
semilogx((1:length(Diff2f))*dt,Diff2f/(fnum),'color',[1 1 1]*.5,'linestyle','-','linewidth',1);
hold off
xlim([.05 30]);
%ylim
ylim([0 25]);
set(gca,'ytick',[0 5 10 15 20]);
%set(gca,'xticklabel',[ ]); %set(gca,'ytick',[ ]);
set(gca,'xtick',10.^([-1 0 1]));
set(gca,'yaxislocation','right');
set(gca,'xaxislocation','top');
%grid on
set(gca,'xticklabel',[]);
xlabel('Error Growth','fontsize',10,'interpreter','tex');
%xlabel('Time (mtu)','fontsize',10)
%box on
line([.06 .09],[1 1]*.835*max(ylim),'color',[1 1 1]*.5);
text(.10,.835*max(ylim),'$|{\bf x}-{\bf x}_{12\:bit}''|(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');
line([.06 .09],[1 1]*.635*max(ylim),'color',[1 1 1]*0);
text(.10,.635*max(ylim),'$|{\bf x}-{\bf x}''|(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');



axes('position',[bxlf+3*bxwd+3*bxpw bxbt+1*bxht+1*bxpu bxwd bxht])
semilogx((1:length(Diffpf))*dt,Diffpf/fnum,'color','k','linestyle','-','linewidth',1); hold on
semilogx((1:length(Diff2f))*dt,Diff3f/fnum,'color',[1 1 1]*.5,'linestyle','-','linewidth',1);
hold off
xlim([.05 30]);
ylim([0 25]);
set(gca,'ytick',[0 5 10 15 20]);
set(gca,'xticklabel',[ ]); %set(gca,'ytick',[ ]);
set(gca,'xtick',10.^([-1 0 1]));
set(gca,'yaxislocation','right');
line([.06 .09],[1 1]*.835*max(ylim),'color',[1 1 1]*.5);
text(.10,.835*max(ylim),'$|{\bf x}-{\bf x}_{14\:bit}''|(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');

%set(gca,'yaxislocation','right');
%grid on

axes('position',[bxlf+3*bxwd+3*bxpw bxbt+0*bxht+0*bxpu bxwd bxht])
semilogx((1:length(Diffpf))*dt,Diffpf/fnum,'color','k','linestyle','-','linewidth',1); hold on
semilogx((1:length(Diff2f))*dt,Diff4f/fnum,'color',[1 1 1]*.5,'linestyle','-','linewidth',1);
hold off
xlim([.05 30]);
ylim([0 25]);
%set(gca,'xticklabel',[ ]); %set(gca,'ytick',[ ]);
set(gca,'ytick',[0 5 10 15 20]);
set(gca,'xtick',10.^([-1 0 1]));
set(gca,'yaxislocation','right');
%set(gca,'yaxislocation','right');
%grid on
xlabel('Time (mtu)','fontsize',10)
line([.06 .09],[1 1]*.835*max(ylim),'color',[1 1 1]*.5);
text(.10,.835*max(ylim),'$|{\bf x}-{\bf x}_{16\:bit}''|(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');




%draw bit boxes
axes('position',[bxlf-2*bxpw bxbt+2*bxht+2.2*bxpu 1.2*bxpw .98*bxht]);
bxclr = ones(1,Numbits);
bxclr(Numbits-4:Numbits-2)=.4;
bxclr(1:bits2)=.4;
%bxclr(bits2-9)=.7;
%drawbitbox3(Numbits,bxclr);
%set limits
numbits=Numbits;
boxcolors = bxclr;
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
text(-8,numbits+1+.4,'Sign','interpreter','tex','fontsize',6,'horizontalalignment','left');
%text(-2.5,numbits+1+.4,'Sign','interpreter','tex','fontsize',6,'horizontalalignment','center');
text(-8,numbits-6+.4,'Exp.','interpreter','tex','fontsize',6,'horizontalalignment','left');
text(-5,numbits-27+.4,'Fraction','interpreter','tex','fontsize',6,'horizontalalignment','left','rotation',90);
% text(-0.25,(23-(32-numbits))/2+.4,'F','interpreter','latex','fontsize',7,'horizontalalignment','center');
text(-26,numbits-12,'12 bit','interpreter','tex','fontsize',10,'fontweight','normal','horizontalalignment','left','rotation',0);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
axis off


%draw bit boxes
axes('position',[bxlf-2*bxpw bxbt+bxht+1.3*bxpu 1.2*bxpw .98*bxht]);
bxclr = ones(1,Numbits);
bxclr(Numbits-4:Numbits-2)=.4;
%bxclr(1:bits3-10)=.7;
bxclr(1:bits3)=.4;
numbits=Numbits;
boxcolors = bxclr;
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
text(-8,numbits+1+.4,'Sign','interpreter','tex','fontsize',6,'horizontalalignment','left');
%text(-2.5,numbits+1+.4,'Sign','interpreter','tex','fontsize',6,'horizontalalignment','center');
text(-8,numbits-6+.4,'Exp.','interpreter','tex','fontsize',6,'horizontalalignment','left');
text(-5,numbits-27+.4,'Fraction','interpreter','tex','fontsize',6,'horizontalalignment','left','rotation',90);
text(-26,numbits-12,'14 bit','interpreter','tex','fontsize',10,'fontweight','normal','horizontalalignment','left','rotation',0);
% text(-0.25,(23-(32-numbits))/2+.4,'F','interpreter','latex','fontsize',7,'horizontalalignment','center');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
axis off


%draw bit boxes
axes('position',[bxlf-2*bxpw 1.01*bxbt 1.2*bxpw .98*bxht]);
bxclr = ones(1,Numbits);
bxclr(Numbits-4:Numbits-2)=.4;
bxclr(1:bits4)=.4;
%bxclr(bits4-9)=.7;
numbits=Numbits;
boxcolors = bxclr;
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
text(-8,numbits+1+.4,'Sign','interpreter','tex','fontsize',6,'horizontalalignment','left');
%text(-2.5,numbits+1+.4,'Sign','interpreter','tex','fontsize',6,'horizontalalignment','center');
text(-8,numbits-6+.4,'Exp.','interpreter','tex','fontsize',6,'horizontalalignment','left');
text(-5,numbits-27+.4,'Fraction','interpreter','tex','fontsize',6,'horizontalalignment','left','rotation',90);
% text(-0.25,(23-(32-numbits))/2+.4,'F','interpreter','latex','fontsize',7,'horizontalalignment','center');
text(-26,numbits-12,'16 bit','interpreter','tex','fontsize',10,'fontweight','normal','horizontalalignment','left','rotation',0);
%text(-16,numbits-16,'16-bit','interpreter','tex','fontsize',10,'fontweight','normal','horizontalalignment','left','rotation',90);

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
axis off





set(gcf,'papersize',[16.5 20.32/2]);
set(gcf,'paperposition',[0 0 16.5 20.32/2]);

saveas(gcf,[savedir 'l63fig4.pdf']);




%% Fig 3

set(0,'defaultaxesfontsize',9);

Numbits = 22;

figure(3)
if dockfigs
    set(gcf,'windowstyle','docked');
else
    set(gcf,'visible','off');
end

clf
set(gcf,'color','w');

%draw bit boxes
axes('position',[.005 .085 .05 .87]);
drawbitbox2(Numbits,0);


axis off
axes('position',[.06 .085 .7 .87]);

clc


timedt = (timeTN(2:end) - timeTN(1:end-1))';
Tinfoy = Tinfox;
Tinfoy(Tinfoy<5e-5)=0;
suminfo = Tinfoy(:,1:end-1)*timedt;
%plot(suminfo,'*');

tmp = real(log(Tinfox(:,1:end-5)));
tmp2 = flipud(tmp);
[tn,tm] = size(tmp2);
tmp3 = ones(tn+2,tm)*-inf;
tmp3(1,:) = tmp2(1,:);
tmp3(3:10,:) = tmp2(2:9,:);
tmp3(12:Numbits+2,:) = tmp2(10:Numbits,:);

tmp=tmp3(1:Numbits+2,:);

tmp(12,15)=1.2*tmp(12,15);
tmp(13,15)=1.2*tmp(13,15);
tmp(14,15)=1.2*tmp(14,15);
tmp(15,15)=1.3*tmp(15,15);
tmp(16,15)=1.4*tmp(16,15);


h = imagesc(tmp); 
cb = colorbar;
pos = get(cb,'position');
pos(1) = 1.32*pos(1);
pos(3) = .7*pos(3);
set(cb,'position',pos);

caxis([-11 1]);
mygr = .05;
cm = zeros(64:3);
cm (:,1) = linspace(mygr,1,64);
cm (:,2) = cm (:,1);
cm (:,3) = cm (:,1);
colormap(flipud(cm));

set(gca,'xtick',[])
set(gca,'ytick',[]);
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
%axis off

hold on

%draw box lines
line(xlim,[min(ylim) min(ylim)],'color',[1 1 1]*.01);
for ln=1:Numbits+1
     line(xlim,[ln ln]+.5,'color',[1 1 1]*0);
end
line(xlim,[max(ylim) max(ylim)],'color',[1 1 1]*.01);
line([min(xlim) min(xlim)],[0 1]+.5,'color','k');
line([max(xlim) max(xlim)],[0 1]+.5,'color','k');
line([min(xlim) min(xlim)],[2 10]+.5,'color','k');
line([max(xlim) max(xlim)],[2 10]+.5,'color','k');
line([min(xlim) min(xlim)],[11 Numbits+2]+.5,'color','k');
line([max(xlim) max(xlim)],[11 Numbits+2]+.5,'color','k');


%text total info
tbn = 1;
jstr = [' ' num2str(suminfo(24-tbn-1),'%1.1e')];
text(1.15*max(xlim),tbn+.1,jstr,'interpreter','tex','fontsize',10,'horizontalalignment','right');
for tbn=3:10
    if suminfo(25-tbn-1)>0
        jstr = [' ' num2str(suminfo(25-tbn-1),'%1.1e')];
    else
        jstr = ['0      '];
    end
    text(1.15*max(xlim),tbn+.1,jstr,'interpreter','tex','fontsize',10,'horizontalalignment','right');
end
for tbn=12:Numbits+1
    if suminfo(26-tbn)>0
        jstr = [' ' num2str(suminfo(26-tbn),'%1.1e')];
    else
        jstr = ['0      '];
    end
    text(1.15*max(xlim),tbn+.1,jstr,'interpreter','tex','fontsize',10,'horizontalalignment','right');
end



%xlabel
text(1.07*mean(xlim),max(ylim)+.08*range(ylim),'Forecast time $\tau$ (mtu)','interpreter','latex','fontsize',10,'horizontalalignment','center');
%title
text(mean(xlim),min(ylim)-.025*range(ylim),'Bitwise information content $I_b$ decaying with forecast time $\tau$','interpreter','latex','fontsize',10,'horizontalalignment','center');
text(52,.8*mean(ylim),...
    'Log[ $I_{b}(\tau)$ ]','interpreter','latex','fontsize',10,...
    'rotation',-90);
text(max(xlim)+.075*range(xlim),min(ylim)-.025*range(ylim),'$J_b$','interpreter','latex','fontsize',10,'horizontalalignment','center');


axis off
axes('position',[.06 .086 .69 .03]);
semilogx(0,'.w');
xlim([0.05 100]);
set(gca,'ytick',[]);
box off
set(gca,'TickLabelInterpreter', 'tex');

%title('h2');

set(gcf,'papersize',[16.5 .61*20.32]);
set(gcf,'paperposition',[0 0 16.5 .61*20.32]);

saveas(gcf,[savedir 'l63fig3.pdf']);


%% FIG 2


figure(2)
if dockfigs
    set(gcf,'windowstyle','docked');
else
    set(gcf,'visible','off');
end


Numbits=32;

clf

bots = linspace(.03,.79,5);
height = (bots(2)-bots(1))*.9;

lfts = linspace(.1,.7,3)+.02;
width = (lfts(2)-lfts(1))*.925;


axes('position',[lfts(1) bots(5) width height])
bitinfoplot(p1baseline,p2b32_t1,p3b32_t1,q2b32_t1,q3b32_t1,2.5,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]); %ylim([0 .05]);
text(double(mean(xlim)),double(1.1*max(ylim)),'\textbf{Instantaneous} ($\tau=0$)','interpreter','latex','fontsize',10,'horizontalalignment','center');
%text(double(mean(xlim)),double(1.08*max(ylim)),...
%    '{Instantaneous}{(\tau}=0)','interpreter','tex','fontsize',10,'horizontalalignment','center');

axes('position',[lfts(2) bots(5) width height])
bitinfoplot(p1baseline,p2b32_t2,p3b32_t2,q2b32_t2,q3b32_t2,2.5,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);%ylim([0 .05]);
text(double(mean(xlim)),double(1.1*max(ylim)),'\textbf{Short forecast} ($\tau=0.2$)','interpreter','latex','fontsize',10,'horizontalalignment','center');

axes('position',[lfts(3) bots(5) width height])
bitinfoplot(p1baseline,p2b32_t3,p3b32_t3,q2b32_t3,q3b32_t3,2.5,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);%ylim([0 .05]);
text(double(mean(xlim)),double(1.1*max(ylim)),'\textbf{Long forecast} ($\tau=1$)','interpreter','latex','fontsize',10,'horizontalalignment','center');

%legend
line([-.9 -.75]*max(xlim), [1 1]*.87*max(ylim) ,'color','k','linestyle',':')
text(double(-.68*max(xlim)),double(.87*max(ylim)),'$p(x)$', ... 
     'interpreter','latex','horizontalalignment','left', ...
     'fontsize',10)
line([-.9 -.75]*max(xlim), [1 1]*.73*max(ylim) ,'color','k')
text(double(-.68*max(xlim)),double(.73*max(ylim)),'$p(x \:| \:x_b(t-\tau) = 1 )$', ... 
     'interpreter','latex','horizontalalignment','left', ...
     'fontsize',10)
line([-.9 -.75]*max(xlim), [1 1]*.59*max(ylim) ,'color',[1 1 1]*.5)
text(double(-.68*max(xlim)),double(.59*max(ylim)),'$p(x \:| \:x_b(t-\tau) = 0 )$', ... 
     'interpreter','latex','horizontalalignment','left', ...
     'fontsize',10)
% line([-.95 -.77]*max(xlim), [1 1]*.72*max(ylim) ,'color',[1 1 1]*.7)
% text(-.75*max(xlim),.72*max(ylim),'$p(x \:| \:x_0(t) = 0)$', ... 
%     'interpreter','latex','horizontalalignment','left', ...
%     'fontsize',8)
% text(.9*max(xlim),.62*max(ylim),'$I=1.00$', ... 
%     'interpreter','latex','horizontalalignment','right', ...
%     'fontsize',8)

axes('position',[lfts(1) bots(4) width height])
bitinfoplot(p1baseline,p2b31_t1,p3b31_t1,q2b31_t1,q3b31_t1,10,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);

axes('position',[lfts(2) bots(4) width height])
bitinfoplot(p1baseline,p2b31_t2,p3b31_t2,q2b31_t2,q3b31_t2,3,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);

axes('position',[lfts(3) bots(4) width height])
bitinfoplot(p1baseline,p2b31_t3,p3b31_t3,q2b31_t3,q3b31_t3,2,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);

axes('position',[lfts(1) bots(3) width height])
bitinfoplot(p1baseline,p2b24_t1,p3b24_t1,q2b24_t1,q3b24_t1,3,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);

axes('position',[lfts(2) bots(3) width height])
bitinfoplot(p1baseline,p2b24_t2,p3b24_t2,q2b24_t2,q3b24_t2,1.5,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);

axes('position',[lfts(3) bots(3) width height])
bitinfoplot(p1baseline,p2b24_t3,p3b24_t3,q2b24_t3,q3b24_t3,2,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);

axes('position',[lfts(1) bots(2) width height])
bitinfoplot(p1baseline,p2b23_t1,p3b23_t1,q2b23_t1,q3b23_t1,3,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);

axes('position',[lfts(2) bots(2) width height])
bitinfoplot(p1baseline,p2b23_t2,p3b23_t2,q2b23_t2,q3b23_t2,1.5,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);

axes('position',[lfts(3) bots(2) width height])
bitinfoplot(p1baseline,p2b23_t3,p3b23_t3,q2b23_t3,q3b23_t3,2,Bins);
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[]);

axes('position',[lfts(1) bots(1) width height])
bitinfoplot(p1baseline,p2b20_t1,p3b20_t1,q2b20_t1,q3b20_t1,1.5,Bins);
box on; set(gca,'xtick',[-15 0 15]); %set(gca,'xticklabel',[]);

axes('position',[lfts(2) bots(1) width height])
bitinfoplot(p1baseline,p2b20_t2,p3b20_t2,q2b20_t2,q3b20_t2,1.5,Bins);
box on; set(gca,'xtick',[-15 0 15]); %set(gca,'xticklabel',[]);

axes('position',[lfts(3) bots(1) width height])
bitinfoplot(p1baseline,p2b20_t3,p3b20_t3,q2b20_t3,q3b20_t3,2,Bins);
box on; set(gca,'xtick',[-15 0 15]); %set(gca,'xticklabel',[]);


% % draw bit boxes
% bxlf = lfts(1)-.115*width;
% bxwd = .075*width;
axes('position',[bxlf bots(5) bxwd height]);
text(-.3,19,'\textbf{Sign}','interpreter','latex','fontsize',10,'horizontalalignment','center');
text(-.3,15.0,'\textbf{bit}','interpreter','latex','fontsize',10,'horizontalalignment','center');
xlim([0 4]);ylim([0 numbits+2]);axis off

% drawbitbox8(Numbits,32-(32-Numbits));
axes('position',[bxlf bots(4) bxwd height]);
text(-.3,24,'\textbf{1\textsuperscript{st}}','interpreter','latex','fontsize',10,'horizontalalignment','center');
text(-.3,19.5,'\textbf{Exponent}','interpreter','latex','fontsize',10,'horizontalalignment','center');
text(-.3,15.0,'\textbf{bit}','interpreter','latex','fontsize',10,'horizontalalignment','center');
xlim([0 4]);ylim([0 numbits+2]);axis off
% drawbitbox8(Numbits,31-(32-Numbits));
axes('position',[bxlf bots(3) bxwd height]);
text(-.3,24,'\textbf{7\textsuperscript{th}}','interpreter','latex','fontsize',10,'horizontalalignment','center');
text(-.3,19.5,'\textbf{Exponent}','interpreter','latex','fontsize',10,'horizontalalignment','center');
text(-.3,15.0,'\textbf{bit}','interpreter','latex','fontsize',10,'horizontalalignment','center');
xlim([0 4]);ylim([0 numbits+2]);axis off
% drawbitbox8(Numbits,24-(32-Numbits));
axes('position',[bxlf bots(2) bxwd height]);
text(-.3,24,'\textbf{1\textsuperscript{st}}','interpreter','latex','fontsize',10,'horizontalalignment','center');
text(-.3,19.5,'\textbf{Fraction}','interpreter','latex','fontsize',10,'horizontalalignment','center');
text(-.3,15.0,'\textbf{bit}','interpreter','latex','fontsize',10,'horizontalalignment','center');
xlim([0 4]);ylim([0 numbits+2]);axis off
% drawbitbox8(Numbits,23-(32-Numbits));
axes('position',[bxlf bots(1) bxwd height]);
text(-.3,24,'\textbf{3\textsuperscript{rd}}','interpreter','latex','fontsize',10,'horizontalalignment','center');
text(-.3,19.5,'\textbf{Fraction}','interpreter','latex','fontsize',10,'horizontalalignment','center');
text(-.3,15.0,'\textbf{bit}','interpreter','latex','fontsize',10,'horizontalalignment','center');
xlim([0 4]);ylim([0 numbits+2]);axis off
% %drawbitbox(Numbits,20-(32-Numbits));
% drawbitbox8(Numbits,20-(32-Numbits));


set(gcf,'papersize',[17.5 20.32]);
set(gcf,'paperposition',[0 0 16.5 20.32]);

saveas(gcf,[savedir 'l63fig2.pdf']);










%% FIG 1


set(0,'defaultaxesfontsize',10);
set(0,'defaulttextfontsize',10);

figure(1)
if dockfigs
    set(gcf,'windowstyle','docked');
else
    set(gcf,'visible','off');
end
clf


bxlf = .081;
bxbt = .11;
bxwd = .26;
bxht = .38;
bxpw = .017;
bxpu = .02;


%plot time series
axes('position',[bxlf bxbt+1*bxht+1*bxpu 2*bxwd+bxpu bxht])
hold on
plot(T1colonPln,X1colonPln,'k'); 
plot(T1colonPln,X1pcolonPln,'color',[1 1 1]*.5,'linestyle','-'); hold off;
ylim([-1 1]*30)
set(gca,'xtick',[0:5:25]);
set(gca,'xticklabel',[0 5 10 15 20 25]); 
set(gca,'xaxislocation','top');
set(gca,'ytick',[-20 -10 0 10 20]);
grid on
ylabel('$x(t)$','interpreter','latex');
xlabel('Timeseries');
hold off
box on

line(0+[.7 2],[1 1]*24.5,'color','k');
text(2.25,25,'$x(t)$','interpreter','latex','horizontalalignment','left','fontsize',10);
line(5+[.7 2],[1 1]*24.5,'color',[1 1 1]*.5);
text(7.25,25,'$x''(t)$','interpreter','latex','horizontalalignment','left','fontsize',10,'color','k');



%plot bitmap
axes('position',[bxlf bxbt 2*bxwd+bxpu bxht])
blcks = 100;
bln = round(linspace(1,Plotlength,blcks));
imagesc(flipud(1-Bitimage'));
colormap gray;
hold off
box on
%draw lines
for bk = 1:blcks-1;
    line([bk bk]+.5,ylim,'color',[1 1 1]*.7);
end
for bt = 1:32-1;
    line(xlim,[bt bt]+.5,'color',[1 1 1]*.7);
end
set(gca,'ytick',[]);
set(gca,'xtick',linspace(min(xlim),max(xlim),6)); 
set(gca,'xticklabel',[0 5 10 15 20 25]); 
yl = ylim;
xlabel('Bit-image');
tmp = ylabel('$x_b(n\Delta t)$','interpreter','latex');
pos = get(tmp,'position');
pos(1)=-8.8;
set(tmp,'position',pos);
%draw bit labels
axes('position',[bxlf-bxpw bxbt bxpw bxht])
ylim(yl);
xlim([0 1])
line(.3+.7*xlim,[0 0]+.5,'color',[1 1 1]*0);
line(.3+.7*xlim,[31 31]+.5,'color',[1 1 1]*0);
line(.3+.7*xlim,[23 23]+.5,'color',[1 1 1]*0);
 text(.99,33,'Sign','interpreter','tex', ...
     'horizontalalignment','right', ...
     'rotation',0, ...
     'fontsize',8);
% text(.55,27.7,'Exp','interpreter','tex', ...
%     'horizontalalignment','right', ...
%     'rotation',0, ...
%     'fontsize',8);
text(-.3,30,'Exp','interpreter','tex', ...
    'horizontalalignment','right', ...
    'rotation',90, ...
    'fontsize',8);
text(-.3,16.5,'Fraction','interpreter','tex', ...
    'horizontalalignment','right', ...
    'rotation',90, ...
    'fontsize',8);
axis off


% plot attractor
axes('position',[bxlf+2*bxwd+2*bxpw bxbt+1*bxht+1*bxpu bxwd bxht])
plot(X1colonPln2,X3colonPln2,'k','linewidth',.5,'color',[1 1 1]*0);
xlim([-1 1]*22)
ylim([0 1]*50)
set(gca,'xtick',[-15 0 15]); 
set(gca,'yaxislocation','right');
set(gca,'xaxislocation','top');
text(1.4*max(xlim),mean(ylim),'$z$','interpreter','latex','fontsize',10);
xlabel('Attractor','interpreter','tex');



%plot sign bit info
pdf=X1pdf;
pdf(X1pdf==0) = .000000001;
entr1 = sum(-pdf.*log2(pdf));
entr2 = sum(-Xsign1pdf.*log2(Xsign1pdf));
entr3 = sum(-Xsign0pdf.*log2(Xsign0pdf));
X1signInfo = entr1-Psign1*entr2-Psign0*entr3;
axes('position',[bxlf+2*bxwd+2*bxpw bxbt bxwd bxht])
hold on
plot(Bins,Xsign1pdf,'color',[1 1 1]*0,'linestyle','-'); 
plot(Bins,Xsign0pdf,'color',[1 1 1]*.5,'linestyle','-'); 
plot(Bins,X1pdf,'color',[1 1 1]*0,'linestyle',':'); 
box on; set(gca,'xtick',[-15 0 15]); set(gca,'xticklabel',[-15 0 15]);
ylim(double([0 1.1*max([X1pdf Xsign1pdf Xsign0pdf])]));
xlim([-1 1]*22)
hold off
set(gca,'xtick',[]);
yl = double(ylim);
yl(2) = 1.5*yl(2);
ylim(yl);
set(gca,'yaxislocation','right');
set(gca,'ytick',[0 .01 .02 .03 .04]);
set(gca,'xtick',[-15 0 15]);
set(gca,'xticklabel',{'-15' ' 0' '15'});
%legend
line([-.95 -.77]*max(xlim), [1 1]*.92*max(ylim) ,'color','k','linestyle',':')
text(-.75*max(xlim),.92*max(ylim),'$p(x)$', ... 
    'interpreter','latex','horizontalalignment','left', ...
    'fontsize',8)
line([-.95 -.77]*max(xlim), [1 1]*.82*max(ylim) ,'color','k')
text(-.75*max(xlim),.82*max(ylim),'$p(x \:| \:x_1(0) = 1 )$', ... 
    'interpreter','latex','horizontalalignment','left', ...
    'fontsize',8)
line([-.95 -.77]*max(xlim), [1 1]*.72*max(ylim) ,'color',[1 1 1]*.5)
text(-.75*max(xlim),.72*max(ylim),'$p(x \:| \:x_1(0) = 0)$', ... 
    'interpreter','latex','horizontalalignment','left', ...
    'fontsize',8)
text(.9*max(xlim),.92*max(ylim),'$H_x=6.65$', ... 
    'interpreter','latex','horizontalalignment','right', ...
    'fontsize',8)
text(.9*max(xlim),.82*max(ylim),'$H_1=5.65$', ... 
    'interpreter','latex','horizontalalignment','right', ...
    'fontsize',8)
text(.9*max(xlim),.72*max(ylim),'$H_0=5.65$', ... 
    'interpreter','latex','horizontalalignment','right', ...
    'fontsize',8)
text(.9*max(xlim),.62*max(ylim),'$I=1.00$', ... 
    'interpreter','latex','horizontalalignment','right', ...
    'fontsize',8)
xlabel('PDF');
text(1.5*max(xlim),mean(ylim),'$p$','interpreter','latex','fontsize',10);

set(gcf,'papersize',[16.5 20.32/2]);
set(gcf,'paperposition',[0 0 16.5 20.32/2]);



saveas(gcf,[savedir 'l63fig1.pdf']);






%%

% Xi = typecast(X(1,:),'uint32');
% 
% n=1;
% tmp = bitget(Xi(n),1:32)==true;
% 
% clf
% imagesc(tmp);
% title(X(1,n));

%



% figure(2)
% clf
% if dockfigs
%     set(gcf,'windowstyle','docked');
% end
% 
% fdn = 1;
% tv = X(1,fdn:end);
% bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
% bnum=32;
% cv = squeeze(Xb(1:length(tv),bnum));
% bitinfoplot(tv,cv,bins);


%%

%bwid=.3;







%%



%% re run with trimmed bits

%x0 = X2(:,end-2);

%load x0final.mat
load x0f2.mat
x0=x;

%%
%K=single(3);
%a = single([10, 28, 8/3]);  %sigma, rho, beta
%a = [20, 600, 3];  %sigma, rho, beta


load l63longrun.mat X;


%%


load Xsp.mat;

spbox = zeros(100,100,100,'int16');

for i=1:1e7
   
    xind = Xsp(1,i);
    yind = Xsp(2,i);
    zind = Xsp(3,i);
    
    xsp = round(100*(20+xind)/40);
    ysp = round(100*(30+yind)/60);
    zsp = round(100*(0+zind)/60);
    
    spbox(xsp,ysp,zsp)=1;
    
    if mod(i,1000)==0
        i/1e7
    end
    
end

%%

%dt = single(    0.01    );
%steps = single(  100000  );
%psize = single(  0.3     );

K= 3;
a = ([10, 28, 8/3]);  %sigma, rho, beta
%a = [20, 600, 3];  %sigma, rho, beta


dt =             0.01;
steps =          5e7;    %50000000;
psize =          0.3;
TNres =          5     ;
bwid = psize;
pevery =            0.01;

jlen    =       1;

fcasn = round(1000/dt);


%steps = single(  1e5  );

%pln = single(round(30/dt));
%steps = pln;

%allocate memory
% dX = zeros(K,1,'single');
% K1=zeros(K,1,'single');
% K2=zeros(K,1,'single');
% K3=zeros(K,1,'single');
% K4=zeros(K,1,'single');
%X2=zeros(K,steps,'single');
X=zeros(K,fcasn,'single');
Xp=zeros(K,fcasn,'single');
X2p=zeros(K,fcasn,'single');
X3p=zeros(K,fcasn,'single');
X4p=zeros(K,fcasn,'single');


bits2=17;
bits3=15;
bits4=13;

load Xsp.mat;
load x0.mat;

%X2(:,1) = x0f;
X(:,1)  = x;
Xp(:,1) = x + psize*randn(3,1);
X2p(:,1)= Xp(:,1);
X3p(:,1)= Xp(:,1);
X4p(:,1)= Xp(:,1);



Diff0f = zeros(1,fcasn);
Diffpf = zeros(1,fcasn);
Diff2f = zeros(1,fcasn);
Diff3f = zeros(1,fcasn);
Diff4f = zeros(1,fcasn);

tic;
fnum=0;
%run full experiment
nup=single(0);

rspot = zeros(3,1);


tv = Xsp(1,1:end);
%Bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
Bins2 = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
pdfXp = zeros(1,length(Bins));
pdfX2p = zeros(1,length(Bins));
pdfX3p = zeros(1,length(Bins));
pdfX4p = zeros(1,length(Bins));



%pdfXp = pdfXp./sum(pdfXp);


tic;
nup=0;
for j=1:jlen
%for i=1:steps-1
for i=1:fcasn
% 
    Dx = RungeLorenz63(X(:,i),dt,a);
    X(:,i+1)=X(:,i)+Dx;
 
   Dx = RungeLorenz63(Xp(:,i),dt,a);
   Xp(:,i+1)=Xp(:,i)+Dx;
    
   Dx = RungeLorenz63(X3p(:,i),dt,a);
   X3p(:,i+1)=X3p(:,i)+Dx;
   
   Dx = RungeLorenz63(X2p(:,i),dt,a);
   X2p(:,i+1)=X2p(:,i)+Dx;
%    
   Dx = RungeLorenz63(X4p(:,i),dt,a);
   X4p(:,i+1)=X4p(:,i)+Dx;
   
   
   X2p(1,i+1) = trimbits(X2p(1,i+1),bits2);
   X2p(2,i+1) = trimbits(X2p(2,i+1),bits2);
   X2p(3,i+1) = trimbits(X2p(3,i+1),bits2);
%    
   X3p(1,i+1) = trimbits(X3p(1,i+1),bits3);
   X3p(2,i+1) = trimbits(X3p(2,i+1),bits3);
   X3p(3,i+1) = trimbits(X3p(3,i+1),bits3);
   
   X4p(1,i+1) = trimbits(X4p(1,i+1),bits4);
   X4p(2,i+1) = trimbits(X4p(2,i+1),bits4);
   X4p(3,i+1) = trimbits(X4p(3,i+1),bits4);
   
   i/fcasn
   
end

    tv = Xp(1,1:end);
    pdfXp = pdfXp+hist(tv,Bins);
    tv = X2p(1,1:end);
    pdfX2p = pdfX2p+hist(tv,Bins);
    tv = X3p(1,1:end);
    pdfX3p = pdfX3p+hist(tv,Bins);
    tv = X4p(1,1:end);
    pdfX4p = pdfX4p+hist(tv,Bins);
    
    
  %if mod(i,fcasn)==0;
        
        %get forecast difference growth
        fnum = fnum+1;
        %trng = ((i-fcasn+1):i);
        
        trng = 1:fcasn;
                
        Diff0ft = 0;
        Diffpft = 0;
        Diff2ft = 0;
        Diff3ft = 0;
        Diff4ft = 0;
        
        %tmp0=0;
        for nk = 1:3
                         
            Diff0ft = Diff0ft + (X(nk,trng)-X(nk,trng(1))).^2;
            Diffpft = Diffpft + (Xp(nk,trng)-X(nk,trng)).^2;
            Diff2ft = Diff2ft + (X2p(nk,trng)-X(nk,trng)).^2;
            Diff3ft = Diff3ft + (X3p(nk,trng)-X(nk,trng)).^2;
            Diff4ft = Diff4ft + (X4p(nk,trng)-X(nk,trng)).^2;
            
        end

        Diff0f = Diff0f + Diff0ft.^(1/2);
        Diffpf = Diffpf + Diffpft.^(1/2);
        Diff2f = Diff2f + Diff2ft.^(1/2);
        Diff3f = Diff3f + Diff3ft.^(1/2);
        Diff4f = Diff4f + Diff4ft.^(1/2);
                 
        X(:,1) = Xsp(:,randi([1 1e7]))+psize*randn(3,1);
        
        Xp(:,1) = X(:,1) + psize*randn(3,1);
        X2p(:,1)= Xp(:,1);
        X3p(:,1)= Xp(:,1);
        X4p(:,1)= Xp(:,1);

  
    
   if j/jlen > nup
       
       nup = nup+pevery;
        compr = nup/toc;
        estcomp = (1-nup)/compr/60/60;
        disp(['complete: ' num2str(nup) '      finish: ' num2str(estcomp) ' hours']);
       %done = nup/stfLN+(stft-1)/stfLN
   end
   
end


%%


attX4p = X4p;
attX3p = X3p;
attX2p = X2p;
attX = X;
attXp = Xp;


plot(attX4p(1,:),attX4p(3,:),'color',[1 1 1]*.5,'linewidth',.5,'color',[1 1 1]*0);
xlim([-1 1]*22)
ylim([0 1]*50)




%%

pdfXp = pdfXp./sum(pdfXp);
pdfX2p = pdfX2p./sum(pdfX2p);
pdfX3p = pdfX3p./sum(pdfX3p);
pdfX4p = pdfX4p./sum(pdfX4p);


%%

clf
%semilogx((1:length(Diff0f))*dt,Diff0f/fnum,'color','b','linestyle','-','linewidth',1); hold on
semilogx((1:length(Diffpf))*dt,Diffpf/fnum,'color','k','linestyle','-','linewidth',1); hold on
semilogx((1:length(Diff2f))*dt,Diff2f/fnum,'color','g','linestyle','-','linewidth',1);
semilogx((1:length(Diff3f))*dt,Diff3f/fnum,'color','r','linestyle','-','linewidth',1);
semilogx((1:length(Diff4f))*dt,Diff4f/fnum,'color','b','linestyle','-','linewidth',1);

hold off

grid on

xlim([.05 30]);



