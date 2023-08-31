function y_dt=detrend(y,q,i,p,T)
% **This function offers three methods to detrend a time series by removing polynomial time trends.**
%
% **Input Parameters:**
%
% [y] - Original series to be detrended.
%
% [q] - Highest order of time trend to be removed.
%       q=-1 for no-detrending; q=0 for demeaning;
%       q=1 for linear trend detrending
%
% [i] - Detrending method: "1" for one-step detrending;
%       "2" for two-step detrending; '3' for GLS detrending.
%
% [p] - Order of lags in the output detrended time series from lag 0 to p.
%       Default is p=0, i.e., no lag.
%
% [T] - Time value at which the time trend ends. Default is the length of the input series.
%
% **Output Parameters:**
%
% [y_dt] - Detrended time series with polynomial trends removed, including lags
%          from lag 0 to p.
%



%定义默认时间长度和滞后阶数
if nargin==4
    T=size(y,1);
elseif nargin==3
    p=0;
end
T1=size(y,1);
y_dt=zeros(T1-p,p+1);
if q==-1
    y_dt=zeros(T-p,p+1);
    for j=1:p+1
        y_dt(:,j)=y(p+2-j:T-j+1,1);
    end
else
    switch i

        case 1
            D=((T-T1+1+p:T)').^(0:q);
            for j=1:p+1
                y_dt(:,j)=y(p+2-j:T1-j+1,1)-D*pinv(D)*y(p+2-j:T1-j+1,1);
            end
        case 2
            D=((1:T1)').^(0:q);
            theta=pinv(D)*y;
            yy_dt=y-D*theta;
            for j=1:p+1
                y_dt(:,j)=yy_dt(p+2-j:T1-j+1,1);
            end
        case 3
            if q==0
                cbar=-7;
            elseif q==1
                cbar=-15.5;
            end
            D=((1:T1)').^(0:q);
            alpha=1+cbar/T1;
            y_a=y;
            y_a(2:T1)=y(2:end)-alpha*y(1:T1-1);
            D_a=D;
            D_a(2:T1,:)=D(2:end,:)-alpha*D(1:T1-1,:);
            phi=pinv(D_a)*y_a;
            yy_dt=y-D*phi;
            y_dt=zeros(T1-p,p+1);
            for j=1:p+1
                y_dt(:,j)=yy_dt(p+2-j:T1-j+1,1);
            end
    end
end

