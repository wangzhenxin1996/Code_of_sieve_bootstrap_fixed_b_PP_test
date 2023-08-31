function [cov,df]=FixedbHAR(u,type,b,adj)
% **Using non-parametric kernel function method to obtain autocorrelation-heteroskedasticity
% robust (HAR) estimates of variance-covariance matrix for input sequence.**
%
% **Input Parameters:**
%
% [u] - Original sequence for which the variance-covariance matrix needs to be estimated.
%
% [type] - Type of kernel function or non-parametric method used when estimating
%          variance-covariance matrix: "bt" for Bartlett kernel;
%          "pr" for Parzen kernel; "qs" for Quadratic Spectral kernel;
%          "os" for Orthonormal Series method.
%
% [b] - Choice of truncation parameter.
%
% [adj] - Indicates if the result is adjusted: "0" for unadjusted; "1" for adjusted,
%         so that the F-statistic calculated using this estimate follows an F-distribution.
%         Default is unadjusted.
%
% **Output Parameters:**
%
% [cov] - Estimate of variance-covariance matrix.
%
% [df] - Size of the second degree of freedom in the standard F-distribution,
%        adjusted to make the F-statistic using this estimate follow an F-distribution.
%


if nargin==3
    adj=0;
end

T=size(u,2);
p=size(u,1);
cov=zeros(p,p);
if strcmpi(type,'bt')||strcmpi(type,'pr')||strcmpi(type,'qs')
    for i=1:T-1
        k=i/(b*T);
        gamma=u(:,1:T-i)*u(:,1+i:T)';
        if strcmpi(type,'qs')
            a=6*pi*k/5;
            cov=cov+3/a^2*(sin(a)/a-cos(a))*(gamma+gamma');
        elseif k>1
            break
        elseif strcmpi(type,'bt')
            cov=cov+(1-k)*(gamma+gamma');
        else
            if k<=0.5
                cov=cov+(1-6*k^2+6*k^3)*(gamma+gamma');
            else
                cov=cov+2*(1-k)^3*(gamma+gamma');
            end
        end
    end
    cov=cov+u*u';
    cov=cov/T;
    
    if strcmpi(type,'bt')
        c1=1;
        c2=2/3;
    elseif strcmpi(type,'pr')
        c1=3/4;
        c2=0.5393;
    elseif strcmpi(type,'qs')
        c1=1.25;
        c2=1;
    end
    
    adjust=0.5*(1+b*(c1+(p-1)*c2))+0.5*exp(b*(c1+(p-1)*c2));
    if adj==1
        cov=cov*adjust;
    end
    
    Kstar=max([floor(1/(b*c2)) p]);
    if strcmpi(type,'bt')
        df=Kstar;
    elseif strcmpi(type,'pr') || strcmpi(type,'qs')
        df=Kstar-p+1;
    end
    
elseif strcmpi(type,'os')
    f=fft(u,[],2);
    F1=real(f);
    F2=imag(f);
    K=b;
    cov=(F1(:,2:K/2+1)*F1(:,2:K/2+1)'+F2(:,2:K/2+1)*F2(:,2:K/2+1)')*2/(T*K);
    if adj==1
        cov=cov*(K/(K-p+1));
    end
    df=K-p+1;
end
end
