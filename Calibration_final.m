function data=Calibration_final(data)
  
% Data is a scalar structure with three fields: fura340, fura 380 and RCaMP, the column for each field 
% represents its temporal evolution and each row represents a different sample. The output data structure has 
% several new fields which are of interest: halfmaxF halfmaxR maxF maxR utime velocityF velocityR

data.fura340=(data.fura340)';
data.fura380=(data.fura380)';
data.rcamp=(data.rcamp)'; 
n=length(data.fura340(:,1))-4;

% Calibration settings for Fura2AM
beta=12.286;
Rmax=12.811;
Rmin=0.394;
interval=140; % Interval of time before and after peak syncronization


% Close all figures
close all


% Calibration of fura and normalization of RCaMP signal by its baseline

for i=1:n
    ratio(i,:)=data.fura340(i,:)./data.fura380(i,:);
    
    data.sCacyto(i,:)= smooth(0.224*beta*((ratio(i,:)-Rmin)./(Rmax - ratio(i,:))));
    data.snormRCAMP(i,:)= smooth (data.rcamp(i,:)/mean(data.rcamp(i,11:51)));
  
    data.Cacyto(i,:)= (0.224*beta*((ratio(i,:)-Rmin)./(Rmax - ratio(i,:))));
    data.normRCAMP(i,:)= (data.rcamp(i,:)/mean(data.rcamp(i,11:51)));
    
    [a, positionF]=max(data.sCacyto(i,:));
    b=mean(data.sCacyto(i,1:11));
    c=((a-b)/2)+b;
    [syncro, syncropointF]= min(abs(data.Cacyto(i,1:positionF)-c));
    data.halfmaxF(i)=(syncropointF-1)*0.2;
    data.maxF(i)=a;
    data.velocityF(i)= max(gradient(data.sCacyto(i,:))/0.2);
    
    
    [d, positionR]=max(data.snormRCAMP(i,:));
    e=mean(data.snormRCAMP(i,1:11));
    f=((d-e)/2)+e;
    [syncroR, syncropointR]= min(abs(data.normRCAMP(i,1:positionR)-f));
    data.halfmaxR(i)=(syncropointR-1)*0.2;
    data.maxR(i)=d;
    data.velocityR(i)= max(gradient(data.snormRCAMP(i,:))/0.2); 
    
    
    
    %uncoupling time (utime)
    data.utime(i)= data.halfmaxR(i)-data.halfmaxF(i);
  
    %syncro plots
    data.Cacytosyncro(i,:)=(data.Cacyto(i,syncropointF-interval:syncropointF+interval));
    data.normRCAMPsyncro(i,:)=(data.normRCAMP(i,syncropointF-interval:syncropointF+interval));
    data.sCacytosyncro(i,:)= smooth(data.Cacytosyncro(i,:));
    data.snormRCAMPsyncro(i,:)= smooth(data.normRCAMPsyncro(i,:));
   endfor


% Plots of the processed data  
   
t=0.2*[-interval:1:interval];
figure
ax1=plot(t,data.sCacytosyncro);
for i=1:n
  hcol(i,:)=get(ax1(i),'Color');
  endfor
hold on
ax2=plot(t,gradient(data.sCacytosyncro));
for i=1:n
  set(ax2(i),'Color',hcol(i,:));
  endfor
for i=1:n
  ccc(i,:)=['Cell ',num2str(i)]';
  endfor

title('Smooth [Ca]cyto and gradient')
legend(ccc,'Location','Northwest')
figure
ax1=plot(t, data.snormRCAMPsyncro);
for i=1:n
  set(ax1(i),'Color',hcol(i,:));
  endfor
hold on
ax1=plot(t,gradient(data.snormRCAMPsyncro));
for i=1:n
  set(ax1(i),'Color',hcol(i,:));
  endfor
legend(ccc,'Location','northwest')
title(' Rcamp syncro')
figure
[ax,h1,h2]= plotyy(t,data.Cacytosyncro, t, data.normRCAMPsyncro);


set(h1,'Marker','o','MarkerSize',5)
legend(h1,ccc,'Location','Northwest')
set(h2,'Marker','x','MarkerSize',5)
title('[Ca]cyto and Rcamp syncro')
for i=1:n
  set(h2(i),'Color',get(h1(i),'Color'));
  endfor
set(ax(2),'YLim',[0,10],'YTickMode','auto');
set(ax(1),'YTickMode','auto','YTick',0:0.2:1.5);

% mean stuff
for j=1:length(data.Cacytosyncro(1,:))
  mplotF(j)=mean(data.Cacytosyncro(:,j));
  splotF(j)=std(data.Cacytosyncro(:,j))/n^.5;
  mplotR(j)=mean(data.normRCAMPsyncro(:,j));
  splotR(j)=std(data.normRCAMPsyncro(:,j))/n^.5;
  endfor
figure
errorbar(t,mplotF,splotF)
hold on
errorbar(t,mplotR,splotR)
D

end