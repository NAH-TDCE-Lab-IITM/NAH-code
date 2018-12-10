clc;

% sumne maadu
data1=Untitled(:,1:30);

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
%for i=1:1:31000
a=data(19000,:);

a1=reshape(a,[6,5]);

a2=a1';

% sensitivity matrix (v/pa)

sens=[0.124 0.124 0.1126 0.1126 0.1126 0.1126;
      0.124 0.1126 0.124 0.124 0.124 0.124;
      0.1126 0.124 0.1126 0.1126 0.1126 0.1126;
      0.1126 0.1126 0.1126 0.1126 0.1126  0.1126;
      0.1126 0.1126 0.1126 0.1126 0.1126 0.1126];
  
  % in pascal 
  a3=a2./sens;
  
  % in db
  
  a4=20*log10(a3/2*100000);
  
contourf(a2)
%end