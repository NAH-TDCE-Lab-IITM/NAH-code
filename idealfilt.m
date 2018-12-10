function [out,ftout,mask] = idealfilt(im,edges,type);

warning off

if nargin<3
    type = 'lp';
end

im=double(im);
mask=zeros(size(im));
[M,N] = size(im);
rM=round(M/2)+1;
rN=round(N/2)+1;
if length(edges) < 2
    [x,y]=meshgrid(1:size(im,1),1:size(im,2));
    ind=((x-rM).^2+(y-rN).^2 < edges^2)';
    mask(ind)=1;
else
    mask(rM-edges(1):rM+edges(1),rN-edges(2):rN+edges(2))=1;
end
mask=fftshift(mask);

if strcmp(type,'hp')
    mask=~mask;
end

if strcmp(type,'gausslp')
    h=fpecial('gaussian',[M,N],edges);
    mask=abs(ifft2(h));
end

if strcmp(type,'laplacian')
    mask=fftshift(-4*pi^2*((x-rM).^2+(y-rN).^2))';
end

 %figure(10)
 %contourf(im)
 %title('Original image')

ft=fft2(im);

% adding a SNR of 20 db
%ftmax=max(max(abs(ft1)));
%SNR=20;
%id=ones(size(ftmax));
%noise=ftmax/100;
%ft=ft1*(1+noise);



  %figure(11)
  %fftshow(fftshift(ft))
  %title('Original FFT')
% % 
%   figure(6)
%   fftshow(fftshift(mask))
%   title('Mask')

 ftout=mask.*fft2(im); 

% no mask

  %ftout=ft;

  %figure(12)
  %fftshow(ftout)
  %title('Output FFT')

out=real(ifft2(ftout));

 % figure(8)
 % contourf(out);
 % title('Filtered image');
% 
% figure(9)
% fftshow(abs(ifft2(mask)));
% title('Spatial filter corresponding to mask');
%%%%%%


    