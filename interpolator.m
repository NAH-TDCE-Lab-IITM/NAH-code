% interpolator

n=4;
a=[1 2 3 4 5;6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20];

[p,q]=size(a);

l=n*(p-1)+1;
m=n*(q-1)+1;
b=zeros(l,m);

for i=0:1:p-1
    for j=0:1:q-1
        b(1+n*i,1+n*j)=a(i+1,j+1);
    end
end


% interpolation along x axis
for k=1:1:p
for j=1:1:q-1
for i=1:1:n-1
    b(1+n*(k-1),i+1+n*(j-1))=(b(1+n*(k-1),1+n*(j-1))*(n-i)+i*b(1+n*(k-1),n+1+n*(j-1)))/n;
end
end
end

% interpolation along y axis
for k=1:1:q
for j=1:1:p-1
for i=1:1:n-1
    b(i+1+n*(j-1),1+n*(k-1))=(b(1+n*(j-1),1+n*(k-1))*(n-i)+i*b(n+1+n*(j-1),1+n*(k-1)))/n;
end
end
end

% interpolation in between



