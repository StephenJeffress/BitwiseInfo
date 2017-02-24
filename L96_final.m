

clear;   %clf; %close all;
clc

K8 = int32(8);
K16 = int32(16);
K32 = int32(32);
K64 = int32(64);

%dt = double(.0005);
%dti = fi(dt,1,16,10);
%dtd = double(dti);
%dts = single(dti);

dtd = .004;
dtd8 = dtd;
dts16  = dtd8/2;
dti32  = fi(dtd8/4,1,16,10);
dtd32 = double(dtd8/4);
dtd64  = dtd8/8;

psizeu=5;
nbias=3.6;

%integration length and update interval
steps = single(  1e6 );
pevr = .0001;

seed = randi(1000);

rng(seed)


rsize8 = 0;
rsize16 = 0;
rsize32 = 0;
rsize64 = 0;

psize1 = 0;
psize16 = 4;
psize32 = 4;
psize64 = 4;
% 

psize64n = .3;


load y8o.mat;
y16           = zeros(K16,1,'single');
y32           = zeros(K32,1,'single');
y64           = zeros(K64,1,'single');
y16(1:2:end)  = y8;
y16(2:2:end)  = nbias+psizeu*(randn(K16/2,1));
y32(1:2:end) = double(y16);
y32(2:2:end) = nbias+psizeu*(randn(K32/2,1));
y64(1:2:end)  = double(y32);
y64(2:2:end)  = nbias+psizeu*(randn(K64/2,1));
        


%Allocate
F = 20;
T = single(dtd:dtd:(dtd*steps));
Y8 = zeros(K8,steps,'single');
Y16 = zeros(K16,steps,'single');
Y32 = zeros(K32,steps,'single');
Y64 = zeros(K64,steps,'single');
Y64p = zeros(K64,steps,'single');

%for perturbation differences
fcast = 5;
fcasn = ceil(fcast/dtd);
fnum = 0;
Diff64 = zeros(1,fcasn);
Diff32 = zeros(1,fcasn);
Diff16 = zeros(1,fcasn);
Diff8 = zeros(1,fcasn);

%get int16 variables
dto2i = fi(dtd32/2,1,16,10);
one6i = fi(1/6,1,16,10);
two6i = fi(2/6,1,16,10);
one2i = fi(1/2,1,16,10);
F32i = fi(F,1,16,10);
Fdti = fi(F*dtd32,1,16,10);
y32i = fi(y32,1,16,10);

%set floating F to equal int16 F*dt
F8 = F*ones(K8,1,'double');
F16 = F*ones(K16,1,'double');
F32 = F*ones(K32,1,'double');
F64 = F*ones(K64,1,'double');



mul16 = fimath('ProductMode','SpecifyPrecision',...
		'ProductWordLength',16,'ProductFractionLength',10);
    




tic
nup=0;


for tt = 1:steps
    
    
    %K=8 in 64 bit double
    Y8(:,tt)=single(y8);
    k1e = -[y8(2:K8); y8(1)].*([y8(3:K8); y8(1:2)]-[y8(K8);y8(1:K8-1)])-y8+F8;
    y1e = y8+k1e*dtd8/2;
    k2e = -[y1e(2:K8); y1e(1)].*([y1e(3:K8); y1e(1:2)]-[y1e(K8);y1e(1:K8-1)])-y1e+F8;
    y2e = y8+k2e*dtd8/2;
    k3e = -[y2e(2:K8); y2e(1)].*([y2e(3:K8); y2e(1:2)]-[y2e(K8);y2e(1:K8-1)])-y2e+F8;
    y3e = y8+k3e*dtd8;
    k4e = -[y3e(2:K8); y3e(1)].*([y3e(3:K8); y3e(1:2)]-[y3e(K8);y3e(1:K8-1)])-y3e+F8;  
    y8 = y8+dtd8/6*(k1e+2*k2e+2*k3e+k4e);
    
    y8 = y8+rsize8*randn(K8,1);
    
    %K=16 in 32 bit single
    Y16(:,tt)=y16;
    for d2=1:2
        k1s = -[y16(2:K16); y16(1)].*([y16(3:K16); y16(1:2)]-[y16(K16);y16(1:K16-1)])-y16+F16;
        y1s = y16+k1s*dts16/2;
        k2s = -[y1s(2:K16); y1s(1)].*([y1s(3:K16); y1s(1:2)]-[y1s(K16);y1s(1:K16-1)])-y1s+F16;
        y2s = y16+k2s*dts16/2;
        k3s = -[y2s(2:K16); y2s(1)].*([y2s(3:K16); y2s(1:2)]-[y2s(K16);y2s(1:K16-1)])-y2s+F16;
        y3s = y16+k3s*dts16;
        k4s = -[y3s(2:K16); y3s(1)].*([y3s(3:K16); y3s(1:2)]-[y3s(K16);y3s(1:K16-1)])-y3s+F16;  
        y16 = y16+dts16/6*(k1s+2*k2s+2*k3s+k4s);
        y16 = y16+rsize16*randn(K16,1);
    end

    
    %save float version of y32i
    
    Y32(:,tt)=single(y32i);
    for d4 = 1:4
    %step forward using integer arithmetic only
    ydti = mpy(mul16, y32i, dti32);
    k1i = accumpos(mpy(mul16,(-[y32i(2:K32); y32i(1)]),...
           accumpos([ydti(3:K32); ydti(1:2)],-[ydti(K32);ydti(1:K32-1)])), ...
           accumpos(-ydti,Fdti));
    
    y1i = accumpos(y32i,mpy(mul16,k1i,one2i));
    y1dti = mpy(mul16, y1i, dti32);
    
    k2i = accumpos(mpy(mul16,(-[y1i(2:K32); y1i(1)]),...
           accumpos([y1dti(3:K32); y1dti(1:2)],-[y1dti(K32);y1dti(1:K32-1)])), ...
           accumpos(-y1dti,Fdti));
    
    y2i = accumpos(y32i,mpy(mul16,k2i,one2i));
    y2dti = mpy(mul16, y2i, dti32);
    
    k3i = accumpos(mpy(mul16,(-[y2i(2:K32); y2i(1)]),...
           accumpos([y2dti(3:K32); y2dti(1:2)],-[y2dti(K32);y2dti(1:K32-1)])), ...
           accumpos(-y2dti,Fdti));
   
    y3i = accumpos(y32i,k3i);
    y3dti = mpy(mul16, y3i, dti32);
    
    k4i = accumpos(mpy(mul16,(-[y3i(2:K32); y3i(1)]),...
           accumpos([y3dti(3:K32); y3dti(1:2)],-[y3dti(K32);y3dti(1:K32-1)])), ...
           accumpos(-y3dti,Fdti));
    
    y32i = accumpos(y32i,...
           accumpos(mpy(mul16,one6i,k1i),...
           accumpos(mpy(mul16,two6i,k2i),...
           accumpos(mpy(mul16,two6i,k3i),...
                    mpy(mul16,one6i,k4i)))));   
    end
    
%     Y32(:,tt)=single(y32);
%     for d4 = 1:4
%         k1d = -[y32(2:K32); y32(1)].*([y32(3:K32); y32(1:2)]-[y32(K32);y32(1:K32-1)])-y32+F32;
%         y1d = y32+k1d*dtd32/2;
%         k2d = -[y1d(2:K32); y1d(1)].*([y1d(3:K32); y1d(1:2)]-[y1d(K32);y1d(1:K32-1)])-y1d+F32;
%         y2d = y32+k2d*dtd32/2;
%         k3d = -[y2d(2:K32); y2d(1)].*([y2d(3:K32); y2d(1:2)]-[y2d(K32);y2d(1:K32-1)])-y2d+F32;
%         y3d = y32+k3d*dtd32;
%         k4d = -[y3d(2:K32); y3d(1)].*([y3d(3:K32); y3d(1:2)]-[y3d(K32);y3d(1:K32-1)])-y3d+F32;  
%         y32 = y32+dtd32/6*(k1d+2*k2d+2*k3d+k4d);
%         y32 = y32+rsize32*randn(K32,1);
%     end
    
    
    %k=64, 64 bit double, truth
    Y64(:,tt)=single(y64);
    for d8 = 1:8
        k1d = -[y64(2:K64); y64(1)].*([y64(3:K64); y64(1:2)]-[y64(K64);y64(1:K64-1)])-y64+F64;
        y1d = y64+k1d*dtd64/2;
        k2d = -[y1d(2:K64); y1d(1)].*([y1d(3:K64); y1d(1:2)]-[y1d(K64);y1d(1:K64-1)])-y1d+F64;
        y2d = y64+k2d*dtd64/2;
        k3d = -[y2d(2:K64); y2d(1)].*([y2d(3:K64); y2d(1:2)]-[y2d(K64);y2d(1:K64-1)])-y2d+F64;
        y3d = y64+k3d*dtd64;
        k4d = -[y3d(2:K64); y3d(1)].*([y3d(3:K64); y3d(1:2)]-[y3d(K64);y3d(1:K64-1)])-y3d+F64;  
        y64 = y64+dtd64/6*(k1d+2*k2d+2*k3d+k4d);
        %y64 = y64+rsize64*randn(K64,1);
    end
    
    %y64           = y64+psize64n*(randn(K64,1));
    
%     Y64p(:,tt)=single(y64p);
%     for d8 = 1:8
%         k1d = -[y64p(2:K64); y64p(1)].*([y64p(3:K64); y64p(1:2)]-[y64p(K64);y64p(1:K64-1)])-y64p+F64;
%         y1d = y64p+k1d*dtd64/2;
%         k2d = -[y1d(2:K64); y1d(1)].*([y1d(3:K64); y1d(1:2)]-[y1d(K64);y1d(1:K64-1)])-y1d+F64;
%         y2d = y64p+k2d*dtd64/2;
%         k3d = -[y2d(2:K64); y2d(1)].*([y2d(3:K64); y2d(1:2)]-[y2d(K64);y2d(1:K64-1)])-y2d+F64;
%         y3d = y64p+k3d*dtd64;
%         k4d = -[y3d(2:K64); y3d(1)].*([y3d(3:K64); y3d(1:2)]-[y3d(K64);y3d(1:K64-1)])-y3d+F64;  
%         y64p = y64p+dtd64/6*(k1d+2*k2d+2*k3d+k4d);
%     end
    
  
    if mod(tt,fcasn)==0;
        
        %get forecast difference growth
        fnum = fnum+1;
        trng = ((tt-fcasn+1):tt);
       
        for pkn=0:7
            Diff8 = Diff8 + abs(Y8(pkn+1,trng)-Y64(8*pkn+1,trng));
            Diff16 = Diff16 + abs(Y16(2*pkn+1,trng)-Y64(8*pkn+1,trng));
            Diff32 = Diff32 + abs(Y32(4*pkn+1,trng)-Y64(8*pkn+1,trng));
            %Diff64 = Diff64 + abs(Y64p(8*pkn+1,trng)-Y64(pkn+1,trng));
        end



        y16(1:2:end)  = y8;
        y16(2:2:end)  = nbias+psizeu*(randn(K16/2,1));
        
        y32(1:2:end) = double(y16);
        y32(2:2:end) = nbias+psizeu*(randn(K32/2,1));
               
        y64(1:2:end)  = double(y32);
        y64(2:2:end)  = nbias+psizeu*(randn(K64/2,1));
        
        y32i = fi(y32,1,16,10);
      
    end
        
    
    if tt/steps>nup
        nup=nup+pevr;
        compr = nup/toc;
        estcomp = (1-nup)/compr/60/60;
        disp(['complete: ' num2str(nup) '      finish: ' num2str(estcomp) ' hours']);
        
    end
    
end

disp('model done');




%plot time series
figure(1)
clf
pause(.1)
plt = 1;
pln = ceil(plt/dtd);
pk = 1;
plot(T(1:pln),Y8(pk,1:pln),'b'); hold on
plot(T(1:pln),Y16(pk,1:pln),'g'); hold on
plot(T(1:pln),Y32(pk,1:pln),'r'); hold on
plot(T(1:pln),Y64(pk,1:pln),'k'); hold off





%plot diff growth
%close(4)


%%


%%

save l96plotvars.mat dft T Diff8 Diff16 Diff32 fnum

%%

load l96plotvars.mat

%%


dft = 14;
set(0,'defaultaxesfontsize',dft);
set(0,'defaulttextfontsize',dft);




figure(2)
%set(gcf,'windowstyle','docked');
set(gcf,'color','white');
clf
subplot('position',[.1 .15 .86 .73])

% lgrng = round(10.^(linspace(0,log10(length(Diff8)),100)));
% lgrng = lgrng(lgrng>1);
% semilogx(T(lgrng(1:4:end)),Diff8(lgrng(1:4:end))/fnum/8,'.','markersize',25,'color','k'); hold on;
% semilogx(T(lgrng(1:2:end)),Diff16(lgrng(1:2:end))/fnum/8,'.','markersize',20,'color','k'); hold on;
% semilogx(T(lgrng),Diff32(lgrng)/fnum/8,'.','markersize',15,'color','k'); hold on;


semilogx(T(1:length(Diff8)),Diff8/fnum/8,'color','k','linestyle','--','linewidth',1); hold on;
semilogx(T(1:length(Diff16)),Diff16/fnum/8,'color','k','linestyle',':','linewidth',1); 
semilogx(T(1:length(Diff32)),Diff32/fnum/8,'color','k','linestyle','-','linewidth',1); 

plot(.1,5.2079,'.k', 'MarkerSize',20);
plot(.1,3.7199,'.k', 'MarkerSize',20);

text(.098,5.4,'5.21','horizontalalignment','right');
text(.11,3.5,'3.72','horizontalalignment','left');


hold off
%grid on

ylim([-1 10]);
xlim([.004 1]);
grid on
set(gca,'xtick',[.01 .1 1]);

%set(gca,'xticklabel',[.01 .1 1]);
title('Forecast accuracy','fontweight','bold','interpreter','latex');
xlabel('Forecast time (mtu)','interpreter','latex')
ylabel('Forecast error','interpreter','latex')
h = legend('Coarse Exact', 'Medium','Fine Inexact','location','northwest');

set(h,'fontsize',dft);
set(h,'interpreter','latex');


%seed

%%


set(gcf,'papersize',[16.5 20.32/2]);
set(gcf,'paperposition',[0 0 16.5 20.32/2]);

saveas(gcf,'/home/jeffress/Documents/Publications/BitwiseInfo/10Jan2017/forecastacc.pdf');

%%


%%
% plot PDF
figure(3)
clf
bwid = .1;
tv = Y8(:);
bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
pdf = hist(tv,bins);
pdf = pdf./sum(pdf);
plot(bins,pdf,'b');
hold on
% tv = Y16(1,:);
% bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
% pdf = hist(tv,bins);
% pdf = pdf./sum(pdf);
% plot(bins,pdf,'g');
% hold on
% tv = Y32(1,:);
% bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
% pdf = hist(tv,bins);
% pdf = pdf./sum(pdf);
% plot(bins,pdf,'r');
tv = Y64(:);
bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
pdf = hist(tv,bins);
pdf = pdf./sum(pdf);
plot(bins,pdf,'r');
hold on

% tv = Y64p(1,:);
% bins = [-fliplr(bwid/2:bwid:abs(min(tv))+bwid/2) (bwid/2:bwid:abs(max(tv))+bwid/2)];
% pdf = hist(tv,bins);
% pdf = pdf./sum(pdf);
% plot(bins,pdf,'r');
% hold on

hold off



%%

clf
plot(1:8:64,y8,'go'); hold on
plot(1:4:64,y16,'b^'); hold on
plot(1:2:64,y32,'rv'); hold on
plot(y64,'k'); hold off

%%


mens = zeros(1,fcasn/2);


cnt = 700;

for nm = 1:cnt
    
    nm/cnt

%nm = 7;

for i=1:fcasn/2
    
    tmp = Y8(:,nm*fcasn+i); 

    mens(i) = mens(i) + mean(tmp);
    
%     plot(tmp);
%     
%     title(num2str(i*dtd));
%     
%     ylim([-1 1]*20);
%     
%     drawnow
%     pause(.01);

    
end

end

%%





