function y=DGP(T,re,c,delta,rho,phi,var)
if all(rho(:)==0,"all") && all(phi(:)==0,"all")
    u=sqrt(var).*randn(T,re);
else
    ARLags=length(rho);
    MALags=length(phi);
    ar=arima(ARLags,0,MALags);
    ar.AR=num2cell(rho);
    ar.MA=num2cell(phi);
    ar.Constant=0;
    ar.Variance=var;
    u=simulate(ar,T,'NumPaths',re);
end
y=zeros(T+1,re);
D=((1:T)').^[0 1];
for i=2:T+1
    y(i,:)=(1+c/T).*y(i-1,:)+u(i-1,:);
end
y(1,:)=[];
Dt=kron(ones(1,re),D*delta);
y=Dt+y;
end