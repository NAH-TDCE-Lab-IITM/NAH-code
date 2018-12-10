function [y]=DiffHankel(n,k,z)

y=n*besselh(n,k,z)/z-besselh(n+1,k,z);
end
