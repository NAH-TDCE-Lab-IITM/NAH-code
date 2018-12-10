%  NAH at a particular frequency
%  NAH TOOLBOX V1
%  Aero-Acoustics group
%  Thermodynamics and Combustion Laboratory
%  IITMADRAS
%  Author : Chaitanya S K -- email id: chaitanya.acharya.007@gmail.com

clear;
clc;

% calling the data file
data1=dlmread('D:\project\experiment\data\axisymmetric cavity\cavity 1\1.6 npr\1.csv',',',[33 1 31032 30]);

% frequency of the source under study.. vary it according to need
l=9709;

% density of air and speed of sound in air
rho=1.2; c=340;

% Number of mics in axial and circumferential direction to the jet
Maxial=5; Mtheta=6;

% total number of mics 
M=Maxial*Mtheta;

% spacing between the mics in axial and circumferential direction 
daxial=0.1;  % in m
dtheta=2*pi/Mtheta;   % in radian
axial_length=daxial*(Maxial-1);
count=0;   % helping to keep in count with number of reconstructions

actual=235;         % radius in mm, do not change it.. 
reconstruction_radius=150:0.5:160;    % radius in mm, change it as per need

% sampling rate of the microphone data
sampling_rate=31000;

%calculate the means
for i=1:1:M
mean1(i)=mean(data1(:,i));
end

mean1=mean1';

% corrected data
for i=1:1:30
    for j=1:1:31000
 data(j,i)=data1(j,i)-mean1(i);  % it is amplitude
    end
end

%% enter the sensitivity matrix         Units : V/Pa  (rms/rms or amp/amp, no change in the value of sensitivity for harmonic source)

% NOTE : the format is like, mic numbers
%   1   2   3   4   5   6
%   7   8   9   10  11  12  
%   13  14  15  16  17  18  
%   19  20  21  22  23  24  
%   25  26  27  28  29  30
%
%   enter the sensitivity of the mic numbers 31 and 32 seperately
 
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
  
sens=reshape(sens',[M,1]);
sens=sens';

for i=1:1:sampling_rate
    sensdata(i,:)=data(i,:)./sens;
     %sensdata(i,:)=data(i,:);
end

% temporal fourier transform of the data
% have a look at this to decide the frequency of the source to be studied.
for i=1:1:M
   S =(sensdata(:,i));
pw_1(:,i)=fourier(sampling_rate,S);
i
end
pw_10=pw_1;
pw_1=pw_1';

sampling_fourier=sampling_rate/2+1;
pw_2=reshape(pw_1,[Mtheta,Maxial,sampling_fourier]);

% making it 0 to 360 deg i.e., repeating the zero degree at 360 degree
Mtheta1=Mtheta+1;

pw_20=zeros(Mtheta1,Maxial,sampling_fourier);
for i=1:1:sampling_fourier
    for j=1:1:Maxial
        for k=1:1:Mtheta
            pw_20(k,j,i)=pw_2(k,j,i);
        end
        pw_20(Mtheta1,j,i)=pw_2(1,j,i);
    end
end

 for i=1:1:sampling_fourier
     for j=1:1:Mtheta1
         for k=1:1:Maxial
     pw_3(k,j,i)=(pw_20(j,k,i));    % temporal signals, these pressures are amplitudes
         end
     end
 end

%  calculating the rms pressures. Since, each amplitude correspond to a
%  harmonic function p(rms)= p(amp)/sqrt(2)
pw_3=pw_3/sqrt(2);
%  figure(1)
% contourf(abs(pw_3(:,:,l)));
% title('hologram cylinder at r=235mm');
% set(gca,'XTickLabel',char('0','30','60','90','120','150','180','210','240','270','300'))
% set(gca,'YTickLabel',char('0','5','10','15','20','25','30','35','40'))


% interpolation of the frequency domain data

b=3;                              % 2^b-1 many points between the points  
for i=1:1:sampling_fourier
a=pw_3(:,:,i);
pw_4(:,:,i)=interp2(a,b,'cubic');   % cubic interpolation
end

% figure(2)
% contourf(abs(pw_4(:,:,l)));
% title('hologram cylinder at r=235mm');

% zero padding
n_zeropad=15;
pzeropad=padarray(abs(pw_4),[n_zeropad 0]);  % 15 in axial direction front and back each. 0 in theta direction

k1=pzeropad(:,:,l);
% figure(3)
% contourf(abs(pzeropad(:,:,l)));


% 2d fourier transform

im=pzeropad;

% choose the 2d image filter 
% lp - low pass filter
% hp - high pass filter. for more look at the idealfilt script

for i=l:1:l
    [out,ftout,mask]=idealfilt(im(:,:,i),[15 15],'lp');         % control the size of the filter
    %out(:,:,i)=out;
    %ftout(:,:,i)=ftout;
    %mask(:,:,i)=mask;
    i
end

%defining kz
MEFTz=((2^b)*(Maxial-1)+1)+(n_zeropad*2);
MEFTtheta=((2^b)*(Mtheta1-1)+1);       

 kappa=2*pi/c*(1:1:sampling_fourier);
 daxial1=daxial*(Maxial-1)/(MEFTz-1);                         % 0.4 is the total length in axial direction
                ksz=2*pi/(daxial1);         % Sampling frequency in z direction
                kz=(ksz)*linspace(-0.5,0.5,MEFTz);     % defining kz
                
% Evaluating kr
                % kr=zeros(MEFTz,MEFTtheta,sampling_rate/2+1);
                kr=zeros(MEFTz,MEFTtheta);
             %for k=2000:1:2000  
              
                    for i=1:MEFTz
                        kr(i,:)=sqrt(+kappa(l)^2-kz(i)^2);                       
                    end
                    
                    
                    

%%  calculation of parameters at different reconstruction radius                    
for count=1:1:length(reconstruction_radius)            
z=[actual reconstruction_radius(count)];            
% pressure propogator
for i=1:1:MEFTz
    for j=1:1:MEFTtheta
       
Gp(i,j)=besselh(1,1,abs(kr(i,j))*z(2))/besselh(1,1,abs(kr(i,j))*z(1));
    end
end

absGp=abs(Gp);

%pressure propogation
P=abs(Gp).*ftout;
P1=(ifft2(P));
P2=P1(n_zeropad+1:MEFTz-n_zeropad,:);                     % removing zerpad
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

% Pressure in db at one chosen frequency
P3(:,:,count)=abs(P2);
P4(:,:,count)=20*log10(P3(:,:,count)*(10^6)/20);

figure(5)
axes1 = axes('Parent',figure(5),'YTickLabel',{'0','10','20','30','40'},...
    'YTick',[1 9 17 25 33],...
    'XTickLabel',{'0','60','120','180','240','300','360'},...
    'XTick',[1 9 17 25 33 41 49],...
    'Layer','top');
box(axes1,'on');
hold(axes1,'all');
a1=reconstruction_radius(count);
title(sprintf('Pressure reconstructed in dB at r=%0.2f cm',a1));
xlabel('angle in degrees');
ylabel('axial distance in mm');
contour(P4(:,:,count),'LineColor',[0 0 0],'Fill','on','Parent',axes1);
colorbar('peer',axes1);

pause(0.5)
clf
% velocity propogator

% d/dz(Hankel(n,k,z))=n*hankel(z)/z-hankel(n+1,k,z);

% W=kr/(1i*c*k) * diffhankel(kr*r)/hankel(kr*r0)*Pn(r0,kz); 
% Fourier Acoustics page 132
for i=1:1:MEFTz
    for j=1:1:MEFTtheta
       a1=kr(i,j)/(1i*rho*c*kappa(l));
Gv(i,j)=a1*DiffHankel(1,1,abs(kr(i,j)*z(2)))/besselh(1,1,abs(kr(i,j)*z(1)));
    end
end

%velocity propogation
V=abs(Gv).*ftout;
V1=(ifft2(V));
V2=V1(n_zeropad+1:MEFTz-n_zeropad,:);
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
I2=I1(n_zeropad+1:MEFTz-n_zeropad,:);
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

count
end

%% plot in r-z  at particular theta at particular frequency

angle_deg=5; % in degree
angle_radian=angle_deg*pi/180;
plot_theta=round(angle_deg/360*(MEFTtheta-1))+1;
lemon=P4(:,plot_theta,:); % (z,theta,r)
lemon1=size(lemon);
lemon2=reshape(lemon,[lemon1(1), lemon1(3)]);

%contourf(lemon2);
xtick_vector1=linspace(1,length(reconstruction_radius),6);
xTickLabel_vector1a=linspace(min(reconstruction_radius),max(reconstruction_radius),6);
xTickLabel_vector1=num2str((xTickLabel_vector1a)');
figure(10)
axes1 = axes('Parent',figure(10),'YTickLabel',{'0','10','20','30','40'},...
    'YTick',[1 9 17 25 33],...
    'XTickLabel',{xTickLabel_vector1},...
    'XTick',xtick_vector1,...
    'Layer','top');
box(axes1,'on');
hold(axes1,'all');
ylabel('axial distance in mm');
xlabel('radius in mm');
title(sprintf('Pressure reconstructed in dB image at theta=%0.2f degree',angle_deg));
contour(lemon2,'LineColor',[0 0 0],'Fill','on','Parent',axes1);
colorbar('peer',axes1);

%% plot in r-theta at particlar z at particular frequency

z_m=0.35;               % in m     0 to 0.4
MEFTz1=MEFTz-n_zeropad*2;
z_plot=round(z_m*(MEFTz1-1)/axial_length)+1;

waterlemon=P4(z_plot,:,:);
waterlemon1=size(waterlemon);
for i=1:1:waterlemon1(2)
    for j=1:1:waterlemon1(3)
        waterlemon2(i,j)=waterlemon(1,i,j);
    end
end
waterlemon2(MEFTtheta,:)=waterlemon2(1,:);
[r,theta] = meshgrid(reconstruction_radius,linspace(0,2*pi,MEFTtheta));
%surf(r.*cos(theta),r.*sin(theta),waterlemon2');
[x,y,v]=pol2cart(theta,r,waterlemon2);
 figure(11)
%h=contourf(x,y,v);
% Create contour
contour(x,y,v,'LineColor',[0 0 0],'Fill','on');

% Create title
title(sprintf('Pressure reconstructed in dB image at z=%0.2f m',z_m));
colorbar;