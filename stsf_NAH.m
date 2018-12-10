%% stsf

% sumne 2nd time maadu
clear all;
clc;

data1=dlmread('D:\project\experiment\data\stsf_2loudspeaker radial\nah\1000 hz\2.csv',',',[33 1 31032 31]);
%data1=hart(:,1:30);

rho=1.2; c=340;
Maxial=5; Mtheta=6;
daxial=0.1; dtheta=2*pi/Mtheta;

z=[235 105];

sampling_rate=31000;
l=1000;

% take 32 instead of 31, call 32 31 [ only in special cases usage]
  % data2(:,1:30)=data1(:,1:30);
  % lemon1=data1(:,32);
  % data3=[data2,lemon1];

%calculate the means
for i=1:1:31
mean1(i)=mean(data1(:,i));  %[normal cases]
%mean1(i)=mean(data3(:,i));   % when only 32 is needed 
end

mean1=mean1';

% corrected data
for i=1:1:31
    for j=1:1:31000
 data(j,i)=data1(j,i)-mean1(i);             % change here as well for only 32
    end
end

% sensitivity matrix 
 
% sens=[104   204    186    142    163    165;
%       164   177    150    177    107    209;
%       88    170    190    205    148    156;
%       146   156    158    151    200    83;
%       140   126    100    99     156    121]/1000;
  
  
sens=[0.1000  0.1558  0.1436  0.1065  0.1133  0.1100
0.1110  0.1160  0.0880  0.1495  0.1090  0.1770
0.0800  0.1460  0.1590  0.1500  0.1270  0.1180
0.1150  0.1186  0.1105  0.1187  0.1380  0.2170
0.1420  0.1290  0.0950  0.0970  0.1604  0.1117];

%0.1690    0.0910
  
sens=reshape(sens',[30,1]);
sens=sens';
sens(1,31)=0.169;                % check for mic sensitivity here
% sensitivity of mics 31 and 32
%sens(1,31)=0.169;    % mic number 31
%sens(1,31)=0.0809;    % mic number 32

for i=1:1:31000
    sensdata(i,:)=data(i,:)./sens;
end


% since both the reference hava a coherence of 1, only one shold be
% choosen,i.e., 1 mic no 31
% calculating auto-correlation of the inputs G3131

q=sensdata(:,31);
r=sensdata(:,31);
G11=cpsd(q,r,rectwin(31000),0,31000,sampling_rate);
figure(1)
cpsd(q,r,rectwin(31000),0,31000,sampling_rate);

% finding Gpp cross-correlation of the input and the outputs

    for j=1:1:30
        q1=sensdata(:,31);
        r1=sensdata(:,j);
        G(1,j,:)=cpsd(q1,r1,rectwin(31000),0,31000,sampling_rate);
    end

% finding the singular value decomposition of Grr @ every frequency.. since
% there is only one input, no need for singular value decomposition..

% finding the pressure field

for i=1:1:15501
    a1=G11(i);
    b=G(1,:,i);
    
    P(:,i)=[(1/a1)*b]'*(a1^0.5);  % careful
    
end

% re-arranging the pressure field

pw_1=reshape(P,[6,5,15501]);

pw_20=zeros(7,5,15501);
for i=1:1:15501
    for j=1:1:5
        for k=1:1:6
            pw_20(k,j,i)=pw_1(k,j,i);
        end
        pw_20(7,j,i)=pw_1(1,j,i);
    end
end

 for i=1:1:15501
     for j=1:1:7
         for k=1:1:5
     pw_3(k,j,i)=(pw_20(j,k,i));    % temporal signals
         end
     end
 end

figure(2)
contourf(abs(pw_3(:,:,l)));
title('hologram cylinder at r=235mm');

% adding 360 degree values 0deg=360deg

for i=1:1:15501
a=pw_3(:,:,i);
pw_4(:,:,i)=interp2(a,3,'cubic');
end

figure(2)
contourf(abs(pw_4(:,:,l)));
title('hologram cylinder at r=235mm');

% zero padding
pzeropad=padarray(abs(pw_4),[15 0]);

k1=pzeropad(:,:,l);
figure(3)
contourf(abs(pzeropad(:,:,l)));


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
P=abs(Gp).*ftout;
P1=(ifft2(P));
P2=P1(16:48,:);
figure(4)
axes1 = axes('Parent',figure(4),'YTickLabel',{'0','10','20','30','40'},...
    'YTick',[1 9 17 25 33],...
    'XTickLabel',{'0','60','120','180','240','300','360'},...
    'XTick',[1 9 17 25 33 41 49],...
    'Layer','top');
box(axes1,'on');
hold(axes1,'all');
title('Pressure reconstructed in Pa image at r=10.5 cm');
contour(abs(P2),'LineColor',[0 0 0],'Fill','on','Parent',axes1);
colorbar('peer',axes1);

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
V=abs(Gv).*ftout;
V1=(ifft2(V));
V2=V1(16:48,:);
figure(6)
axes1 = axes('Parent',figure(6),'YTickLabel',{'0','10','20','30','40'},...
    'YTick',[1 9 17 25 33],...
    'XTickLabel',{'0','60','120','180','240','300','360'},...
    'XTick',[1 9 17 25 33 41 49],...
    'Layer','top');
box(axes1,'on');
hold(axes1,'all');
title('Velocity reconstructed in m/s image at r=10.5 cm');
contour(abs(V2),'LineColor',[0 0 0],'Fill','on','Parent',axes1);
colorbar('peer',axes1);

% velocity in dB
V3=abs(V2);
V4=20*log10(V3*(10^8)/5);

figure(7)
axes1 = axes('Parent',figure(7),'YTickLabel',{'0','10','20','30','40'},...
    'YTick',[1 9 17 25 33],...
    'XTickLabel',{'0','60','120','180','240','300','360'},...
    'XTick',[1 9 17 25 33 41 49],...
    'Layer','top');
box(axes1,'on');
hold(axes1,'all');
title('Velocity reconstructed in dB image at r=10.5 cm');
contour(V4,'LineColor',[0 0 0],'Fill','on','Parent',axes1);
colorbar('peer',axes1);

% intensity calculation

%Power= real(P*conj(V))/2 or dot product of P and V over all points  
%Intensity is power per unit normal area

I1=abs(P1).*abs(V1);
I2=I1(16:48,:);
figure(8)
axes1 = axes('Parent',figure(8),'YTickLabel',{'0','10','20','30','40'},...
    'YTick',[1 9 17 25 33],...
    'XTickLabel',{'0','60','120','180','240','300','360'},...
    'XTick',[1 9 17 25 33 41 49],...
    'Layer','top');
box(axes1,'on');
hold(axes1,'all');
title('Intensity reconstructed image in W/m^2 at r=10.5 cm');
contour(abs(I2),'LineColor',[0 0 0],'Fill','on','Parent',axes1);
colorbar('peer',axes1);            

I3=abs(I2);
I4=10*log10(I3*(10^12));
figure(9)
axes1 = axes('Parent',figure(9),'YTickLabel',{'0','10','20','30','40'},...
    'YTick',[1 9 17 25 33],...
    'XTickLabel',{'0','60','120','180','240','300','360'},...
    'XTick',[1 9 17 25 33 41 49],...
    'Layer','top');
box(axes1,'on');
hold(axes1,'all');
title('Intensity reconstructed image in dB at r=10.5 cm');
contour(I4,'LineColor',[0 0 0],'Fill','on','Parent',axes1);
colorbar('peer',axes1);   
