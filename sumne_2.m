% sumne 2nd time maadu

clear;
clc;

data1=dlmread('D:\project\experiment\data\axisymmetric cavity\cavity 1\1.6 npr\1.csv',',',[33 1 31032 30]);
%data1=hart(:,1:30);

l=9709;  % frequency of source

rho=1.2; c=340;
Maxial=5; Mtheta=6;
daxial=0.1; dtheta=2*pi/Mtheta;

z=[235 200];

sampling_rate=31000;

%calculate the means
for i=1:1:30
mean1(i)=mean(data1(:,i));
end

mean1=mean1';

% corrected data
for i=1:1:30
    for j=1:1:31000
 data(j,i)=data1(j,i)-mean1(i);
    end
end

% sensitivity matrix 
% sens=[0.1126 0.1126 0.124 0.1 0.1126 0.1126;
%       0.1126 0.1126 0.3 0.1126 0.1126 0.1126;
%       0.1126 0.124 0.1126 0.1126 0.1126 0.1126;
%       0.1126 0.1126 0.1126 0.1126 0.1126  0.1126*0.8;
%       0.1126 0.1126 0.1126 0.1126 0.1126 0.124];

% sens of 31 and 32 are 152 and 60 respectively

% sens=[104   204    186    142    163    165;
%       164   177    150    177    107    209;
%       88    170    190    205    148    156;
%       146   156    158    151    200    83;
%       140   126    100    99     156    121]/1000;

sens=[ 0.01045  0.1901  0.1885  0.1427  0.1578  0.1504
  0.1640   0.1584  0.1445  0.1898  0.1553  0.2191
  0.1006   0.2032  0.1873  0.18    0.144   0.1261
  0.1335   0.1287  0.124   0.144   0.1598  0.112
  0.1776   0.1571  0.1166  0.133   0.217   0.159];

%0.1690  0.0910

  
sens=reshape(sens',[30,1]);
sens=sens';

for i=1:1:31000
    sensdata(i,:)=data(i,:)./sens;
     %sensdata(i,:)=data(i,:);
end

% temporal fourier transform of the data

for i=1:1:30
   S =(sensdata(:,i));
pw_1(:,i)=fourier(sampling_rate,S);
i
end
pw_10=pw_1;
pw_1=pw_1';

pw_2=reshape(pw_1,[6,5,15501]);

% making it 0 to 360 deg i.e., repeating the zero degree at 360 degree

pw_20=zeros(7,5,15501);
for i=1:1:15501
    for j=1:1:5
        for k=1:1:6
            pw_20(k,j,i)=pw_2(k,j,i);
        end
        pw_20(7,j,i)=pw_2(1,j,i);
    end
end

 for i=1:1:15501
     for j=1:1:7
         for k=1:1:5
     pw_3(k,j,i)=(pw_20(j,k,i));    % temporal signals
         end
     end
 end

%  figure(1)
% contourf(abs(pw_3(:,:,l)));
% title('hologram cylinder at r=235mm');
% set(gca,'XTickLabel',char('0','30','60','90','120','150','180','210','240','270','300'))
% set(gca,'YTickLabel',char('0','5','10','15','20','25','30','35','40'))
% interpolating of the frequency domain data

for i=1:1:15501
a=pw_3(:,:,i);
pw_4(:,:,i)=interp2(a,3,'cubic');
end
% 
% figure(2)
% contourf(abs(pw_4(:,:,l)));
% title('hologram cylinder at r=235mm');

% zero padding
pzeropad=padarray(abs(pw_4),[15 0]);

k1=pzeropad(:,:,l);
% figure(3)
% contourf(abs(pzeropad(:,:,l)));


% 2d fourier transform

im=pzeropad;

for i=l:1:l
    [out,ftout,mask]=idealfilt(im(:,:,i),[15 15],'lp');
    %out(:,:,i)=out;
    %ftout(:,:,i)=ftout;
    %mask(:,:,i)=mask;
    i
end

%defining kz
 MEFTz=63; MEFTtheta=49;
 kappa=2*pi/c*(1:1:sampling_rate/2+1);
 
                ksz=2*pi/(daxial);         % Sampling frequency in z direction
                kz=(ksz)*linspace(-0.5,0.5,MEFTz);     % defining kz
                
% Evaluating kr
                % kr=zeros(MEFTz,MEFTtheta,sampling_rate/2+1);
                kr=zeros(MEFTz,MEFTtheta);
             %for k=2000:1:2000  
               k=l;
                    for i=1:MEFTz
                        kr(i,:)=sqrt(+kappa(k)^2-kz(i)^2);                       
                    end                
               
                
            % end
             
% pressure propogator
for i=1:1:MEFTz
    for j=1:1:MEFTtheta
       
Gp(i,j)=besselh(1,1,abs(kr(i,j))*z(2))/besselh(1,1,abs(kr(i,j))*z(1));
    end
end

absGp=abs(Gp);

%pressure propogation
P=(Gp).*ftout;
P1=(ifft2(P));
P2=P1(16:48,:);
% figure(4)
% axes1 = axes('Parent',figure(4),'YTickLabel',{'0','10','20','30','40'},...
%     'YTick',[1 9 17 25 33],...
%     'XTickLabel',{'0','60','120','180','240','300','360'},...
%     'XTick',[1 9 17 25 33 41 49],...
%     'Layer','top');
% box(axes1,'on');
% hold(axes1,'all');
% title('Pressure reconstructed in Pa image at r=10.5 cm');
% contour(abs(P2),'LineColor',[0 0 0],'Fill','on','Parent',axes1);
% colorbar('peer',axes1);

% Pressure in db
P3=abs(P2);
P4=20*log10(P3*(10^6)/20);
figure(5)
axes1 = axes('Parent',figure(5),'YTickLabel',{'0','10','20','30','40'},...
    'YTick',[1 9 17 25 33],...
    'XTickLabel',{'0','60','120','180','240','300','360'},...
    'XTick',[1 9 17 25 33 41 49],...
    'Layer','top');
box(axes1,'on');
hold(axes1,'all');
title('Pressure reconstructed in dB image at r=10.5 cm');
contour(P4,'LineColor',[0 0 0],'Fill','on','Parent',axes1);
colorbar('peer',axes1);

% velocity propogator

% d/dz(Hankel(n,k,z))=n*hankel(z)/z-hankel(n+1,k,z);

% W=kr/(1i*c*k) * diffhankel(kr*r)/hankel(kr*r0)*Pn(r0,kz); 
% Fourier Acoustics page 132
for i=1:1:MEFTz
    for j=1:1:MEFTtheta
       a1=kr(i,j)/(1i*rho*c*kappa(k));
Gv(i,j)=a1*DiffHankel(1,1,abs(kr(i,j)*z(2)))/besselh(1,1,abs(kr(i,j)*z(1)));
    end
end

%velocity propogation
V=(Gv).*ftout;
V1=(ifft2(V));
V2=V1(16:48,:);
% figure(6)
% axes1 = axes('Parent',figure(6),'YTickLabel',{'0','10','20','30','40'},...
%     'YTick',[1 9 17 25 33],...
%     'XTickLabel',{'0','60','120','180','240','300','360'},...
%     'XTick',[1 9 17 25 33 41 49],...
%     'Layer','top');
% box(axes1,'on');
% hold(axes1,'all');
% title('Velocity reconstructed in m/s image at r=10.5 cm');
% contour(abs(V2),'LineColor',[0 0 0],'Fill','on','Parent',axes1);
% colorbar('peer',axes1);

% velocity in dB
V3=abs(V2);
V4=20*log10(V3*(10^8)/5);

% figure(7)
% axes1 = axes('Parent',figure(7),'YTickLabel',{'0','10','20','30','40'},...
%     'YTick',[1 9 17 25 33],...
%     'XTickLabel',{'0','60','120','180','240','300','360'},...
%     'XTick',[1 9 17 25 33 41 49],...
%     'Layer','top');
% box(axes1,'on');
% hold(axes1,'all');
% title('Velocity reconstructed in dB image at r=10.5 cm');
% contour(V4,'LineColor',[0 0 0],'Fill','on','Parent',axes1);
% colorbar('peer',axes1);

% intensity calculation

%Power= real(P*conj(V))/2 or dot product of P and V over all points  
%Intensity is power per unit normal area

%I1=P1.*V1;
I1=abs(P1).*abs(V1);
I2=I1(16:48,:);
% figure(8)
% axes1 = axes('Parent',figure(8),'YTickLabel',{'0','10','20','30','40'},...
%     'YTick',[1 9 17 25 33],...
%     'XTickLabel',{'0','60','120','180','240','300','360'},...
%     'XTick',[1 9 17 25 33 41 49],...
%     'Layer','top');
% box(axes1,'on');
% hold(axes1,'all');
% title('Intensity reconstructed image in W/m^2 at r=10.5 cm');
% contour(abs(I2),'LineColor',[0 0 0],'Fill','on','Parent',axes1);
% colorbar('peer',axes1);            

I3=abs(I2);
I4=10*log10(I3*(10^12));
% figure(9)
% axes1 = axes('Parent',figure(9),'YTickLabel',{'0','10','20','30','40'},...
%     'YTick',[1 9 17 25 33],...
%     'XTickLabel',{'0','60','120','180','240','300','360'},...
%     'XTick',[1 9 17 25 33 41 49],...
%     'Layer','top');
% box(axes1,'on');
% hold(axes1,'all');
% title('Intensity reconstructed image in dB at r=10.5 cm');
% contour(I4,'LineColor',[0 0 0],'Fill','on','Parent',axes1);
% colorbar('peer',axes1);   
