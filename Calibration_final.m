AA=DKOR;
AA.fura340=(AA.fura340)';
AA.fura380=(AA.fura380)';
AA.rcamp=(AA.rcamp)'
n=length(AA.fura340(:,1))-4;
beta=12.286;
Rmax=12.811;
Rmin=0.394;

% Close all figures
close all
% Calibration of fura and RCaMP signal (ssmooth=1 to smooth average data before
% analysis and ssmooth=0 to use raw data)

for i=1:n
    ratio(i,:)=AA.fura340(i,:)./AA.fura380(i,:);
    
    AA.sCacyto(i,:)= smooth(0.224*beta*((ratio(i,:)-Rmin)./(Rmax - ratio(i,:))));
    AA.snormRCAMP(i,:)= smooth (AA.rcamp(i,:)/mean(AA.rcamp(i,11:51)));
  
    AA.Cacyto(i,:)= (0.224*beta*((ratio(i,:)-Rmin)./(Rmax - ratio(i,:))));
    AA.normRCAMP(i,:)= (AA.rcamp(i,:)/mean(AA.rcamp(i,11:51)));
    
    end


interval=140;

for i=1:n
    [a, positionF]=max(AA.sCacyto(i,:));
    b=mean(AA.sCacyto(i,1:11));
    c=((a-b)/2)+b;
    [syncro, syncropointF]= min(abs(AA.Cacyto(i,1:positionF)-c));
    AA.halfmaxF(i)=(syncropointF-1)*0.2;
    AA.maxF(i)=a;
    AA.velocityF(i)= max(gradient(AA.sCacyto(i,:))/0.2);
    
    
    [d, positionR]=max(AA.snormRCAMP(i,:));
    e=mean(AA.snormRCAMP(i,1:11));
    f=((d-e)/2)+e;
    [syncroR, syncropointR]= min(abs(AA.normRCAMP(i,1:positionR)-f));
    AA.halfmaxR(i)=(syncropointR-1)*0.2;
    AA.maxR(i)=d;
    AA.velocityR(i)= max(gradient(AA.snormRCAMP(i,:))/0.2); 
    
    
    
    %uncoupling time (utime)
    AA.utime(i)= AA.halfmaxR(i)-AA.halfmaxF(i);
  
    %syncro plots
    AA.Cacytosyncro(i,:)=(AA.Cacyto(i,syncropointF-interval:syncropointF+interval));
    AA.normRCAMPsyncro(i,:)=(AA.normRCAMP(i,syncropointF-interval:syncropointF+interval));
    AA.sCacytosyncro(i,:)= smooth(AA.Cacytosyncro(i,:))
    AA.snormRCAMPsyncro(i,:)= smooth(AA.normRCAMPsyncro(i,:))
   end

DKOR.final=AA;
 
t=0.2*[-interval:1:interval];
figure
ax1=plot(t,AA.sCacytosyncro);
for i=1:n
hcol(i,:)=get(ax1(i),'Color');
end
hold on
ax2=plot(t,gradient(AA.sCacytosyncro));
for i=1:n
set(ax2(i),'Color',hcol(i,:));
end
for i=1:n
    ccc(i,:)=['Celula ',num2str(i)]';
end

title('Smooth [Ca]cyto and gradient')
legend(ccc,'Location','Northwest')
figure
ax1=plot(t, AA.snormRCAMPsyncro);
for i=1:n
set(ax1(i),'Color',hcol(i,:));
end
hold on
ax1=plot(t,gradient(AA.snormRCAMPsyncro));
for i=1:n
set(ax1(i),'Color',hcol(i,:));
end
legend(ccc,'Location','northwest')
title(' Rcamp syncro')
%plot(t,smooth(AA.Cacytosyncro));
%hold on;
%plot(t,gradient(AA.Cacytosyncro));

figure
[ax,h1,h2]= plotyy(t,AA.Cacytosyncro, t, AA.normRCAMPsyncro);
%set(ax(1),'YLim',[0,1]);
set(h1,'Marker','o','MarkerSize',5)
legend(h1,ccc,'Location','Northwest')
set(h2,'Marker','x','MarkerSize',5)
title('[Ca]cyto and Rcamp syncro')
for i=1:n
set(h2(i),'Color',get(h1(i),'Color'));
end
set(ax(2),'YLim',[0,10],'YTickMode','auto');
set(ax(1),'YTickMode','auto','YTick',0:0.2:1.5);
% mean stuff
for j=1:length(AA.Cacytosyncro(1,:))
mplotF(j)=mean(AA.Cacytosyncro(:,j));
splotF(j)=std(AA.Cacytosyncro(:,j))/n^.5;
mplotR(j)=mean(AA.normRCAMPsyncro(:,j));
splotR(j)=std(AA.normRCAMPsyncro(:,j))/n^.5;
end
figure
errorbar(t,mplotF,splotF)
hold on
errorbar(t,mplotR,splotR)
clearvars ax1 ax2 ccc hcol n1 n beta t j mplotR splotR mplotF splotF Rmax Rmin syncroR ssmooth i ratio a b c d e f h1 h2 positionR positionF interval ax  BB syncro syncropointF syncropointR
%clearvars -except AA WT halfmaxF halfmaxR maxF maxR utime velocityF velocityR

