function [zalpha,p_zalpha,zt,p_zt]=ztest(y,dt,q,m,ic,B,b,type)
%%
% **This function performs the PP unit root test on a time series.**
%
% **Input Parameters:**
%
% [y] - Original series for testing.
%
% [dt] - Detrending method:'1' for one-step detrending;
%     '2' for two-step detrending, '3' for GLS detrending.
%
% [q] - Highest power of deterministic linear time trend to be removed.
%       q=-1 for no-detrending; q=0 for demeaning;
%       q=1 for linear trend detrending
%
% [m] - Method for PP test: "1" for traditional PP test,
%      which uses Newey-West for long-term variance estimation;
%      "2" for fixed-b PP unit root test, PP(fb), proposed by Vogelsang and Wagner (2013),
%      which replaces Newey-West long-run variance estimation with fixed-b HAR long-run variance estimation;
%      "4" for obtaining critical values of PP(fb) test using Sieve Bootstrap method, PP^b(fb).
%
% [B] - Number of times bootstrap samples are constructed in the Bootstrap method.
%
% [b] - Truncation parameter for fixed smoothing parameter HAR long-term variance estimation:
%      when "type" is "bt", "pr", or "qs", b represents the truncation parameter;
%      when "type" is "os", it represents the length of the orthogonal sequence used.
%
% [type] - Type of kernel or non-parametric method used when estimating
%          long-term variance-covariance matrix: "bt" for Bartlett kernel;
%          "pr" for Parzen kernel; "qs" for Quadratic Spectral kernel;
%          "os" for Orthonormal Series method.
%
% **Output Parameters:**
%
% [zalpha] - Statistic obtained from the PP test.
%
% [zt] - Statistic obtained from the PP test.
%
% [p_zalpha] - Corresponding p-value for the [zalpha] statistic.
%
% [p_zt] - Corresponding p-value for the [zt] statistic.
%
% **Reference:**
% Vogelsang, T. J. and Wagner, M. (2013). A Fixed-b Perspective on the Phillips-Perron Unit Root Tests.
% Econometric Theory, 29(3):609-628.




%%

if nargin<8,type='bt';end
if nargin<7,b=0.02;end

%{
b_opt=b;
switch type
    case {'bt','pr','qs'}
        
        if b_opt>=0.5
            b=0.5;
        else
            b=b_opt;
        end
        
    case {'os'}
        
        K_opt=cell(1/b_opt);
        if K_opt<=(p+4)
            b=2*floor((p+4)/2);
        elseif K_opt>(p+4) && K_opt<=T
            b=2*floor(K_opt/2);
        else
            b=2*floor(T/2);
        end
        
end
%}



%%
T=size(y,1);
switch m
    case 1
     lags=floor(4*(T/100)^(2/9));
  %     lags=floor(b*T);
        if q==0
            [~,pvalue,stat,~,~]=pptest(y,'model','ARD','test',{'t2','t1'},'lags',lags);
        elseif q==1
            [~,pvalue,stat,~,~]=pptest(y,'model','TS','test',{'t2','t1'},'lags',lags);
        elseif q==-1
            [~,pvalue,stat,~,~]=pptest(y,'model','AR','test',{'t2','t1'},'lags',lags);
        end
        zalpha=stat(1);
        zt=stat(2);
        p_zalpha=pvalue(1);
        p_zt=pvalue(2);
    case 2
        [zalpha,zt]=mztest(y,dt,q,b,type);
        [p_zalpha,p_zt]=CV(q,dt,zalpha,zt,type,b);
        
    case 3
        [zalpha,zt]=mztest(y,dt,q,b,type);
        [wbmzalpha,wbmzt]=wildbootstrap(y,dt,q,B,b,type);
        p_zalpha=sum(wbmzalpha<=zalpha)/B;
        p_zt=sum(wbmzt<=zt)/B;
        
    case 4
        [zalpha,zt]=mztest(y,dt,q,b,type);
        [sbmzalpha,sbmzt]=sievebootstrap(y,dt,q,B,ic,b,type);
        p_zalpha=sum(sbmzalpha<=zalpha)/B;
        p_zt=sum(sbmzt<=zt)/B;
        
        
end
end

function [sbmzalpha,sbmzt]=sievebootstrap(y,dt,q,B,ic,b,type)

T1=size(y,1);
dty=detrend(y,q,dt,0);
n=T1-1;
MaxLags=floor(12*(n/100)^(1/4));
IC=zeros(MaxLags,1);


Ddty=diff(dty);

ye=Ddty(MaxLags+1:end);
switch ic
    case {'aic','maic'}
        CT=2;
    case {'bic','mbic'}
        CT=log(T1-MaxLags);
end

for p=1:MaxLags
    xe=NaN(length(ye),p);
    xe(:,1)=dty(MaxLags+1:end-1);
        for f=1:p
            xe(:,f+1)=Ddty(MaxLags-f+1:end-f);
        end
    alpha_hat=xe\ye;
    u=ye-xe*alpha_hat;
    ic_1=log(u'*u/(T1-MaxLags));
    switch ic
        case {'aic','bic'}
            tau=0;
        case {'maic','mbic'}
            tau=(alpha_hat(1,1))^2*(dty(MaxLags+1:n)'*dty(MaxLags+1:n))/(u'*u/(T1-MaxLags));
    end
    
    IC(p) = ic_1+CT*(tau+p)/(T1-MaxLags);
end

N=min(IC);
[kc,~]=find(IC==N);





switch dt
    case 1
        DLy=detrend(diff(y),q,dt,kc,T1);
        T2=size(DLy,1);
      
    case {2,3,-1}
        Ly=detrend(y,q,dt,kc);
        DLy=diff(Ly);
        T2=size(DLy,1);
end

X=DLy(:,2:end);
Y=DLy(:,1);
%% OLS-method estimates autocorrelated coefficients
beta_hat=pinv(X)*Y;

%% Yule-Waller method estimates autocorrelated coefficients
%  beta_hat=aryule(Ddty,kc);
%  beta_hat=-1.*beta_hat(2:end)';

e=Y-X*beta_hat;
e=e-mean(e);

%bootstrap
sbmzalpha=zeros(B,1);
sbmzt=zeros(B,1);

drop_num=200;
for i=1:B
    index=unidrnd(T2,T1+drop_num,1);
    bre=e(index);
    BU=zeros(T1+drop_num,1);
    for j=kc+1:T1+drop_num
        BU(j,1)=BU(j-1:-1:j-kc,1)'*beta_hat+bre(j);
    end
    
    BU(1:drop_num)=[];
    By=cumsum(BU);
    
    [sbmzalpha(i,1),sbmzt(i,1)]=mztest(By,dt,q,b,type);
end

end

function [wbmzalpha,wbmzt]=wildbootstrap(y,dt,q,B,b,type)
ydt=detrend(y,q,dt,1);
T=size(ydt,1);
y_dt=ydt(:,1);
y_dt_1=ydt(:,2);
alpha_hat=pinv(y_dt_1)*y_dt;
u_hat=y_dt-alpha_hat*y_dt_1;

wbmzalpha=zeros(B,1);
wbmzt=zeros(B,1);

for j=1:B
    w=randn(T,1);
    wbu=u_hat.*w;
    wby=cumsum(wbu);
    [wbmzalpha(j,1),wbmzt(j,1)]=mztest(wby,dt,q,b,type);
end

end

function [mzalpha, mzt]=mztest(y,dt,q,b,type)
ydt=detrend(y,q,dt,1);
T=size(ydt,1);
y_dt=ydt(:,1);
y_dt_1=ydt(:,2);
alpha_hat=pinv(y_dt_1)*y_dt;
u_hat=y_dt-alpha_hat*y_dt_1;
sigmasq=(1/(T-1))*(u_hat'*u_hat);
t_alpha=(alpha_hat-1)/sqrt(sigmasq*(y_dt_1'*y_dt_1)^(-1));
malpha_hat=alpha_hat+0.5*sigmasq/(T^(-1)*(y_dt_1'*y_dt_1));
mu_hat=y_dt-malpha_hat*y_dt_1;
mcov=FixedbHAR(mu_hat',type,b,0);
mzalpha=T*(alpha_hat-1)-0.5*(mcov-sigmasq)*(T^(-2)*(y_dt_1'*y_dt_1))^(-1);
mzt=sqrt(sigmasq/mcov)*t_alpha-0.5*(mcov-sigmasq)*(mcov*T^(-2)*(y_dt_1'*y_dt_1))^(-1/2);
end

% function [zalpha,zt]=standardztest(y,dt,q)
% ydt=detrend(y,q,dt,1);
% T=size(ydt,1);
% y_dt=ydt(:,1);
% y_dt_1=ydt(:,2);
%
% alpha_hat=pinv(y_dt_1)*y_dt;
% u_hat=y_dt-alpha_hat*y_dt_1;
% sigmasq=(1/(T-3))*(u_hat'*u_hat);
% t_alpha=(alpha_hat-1)/sqrt(sigmasq*(y_dt_1'*y_dt_1)^(-1));
% lags=floor(4*(T/100)^(2/9));
% b=(lags+1)/T;
% cov=FixedbHAR(u_hat','bt',b,0);
% sigmasq=(1/(T))*(u_hat'*u_hat);
% zalpha=(T)*(alpha_hat-1)-0.5*(cov-sigmasq)*((T)^2)./(y_dt_1'*y_dt_1);
% zt=sqrt(sigmasq/cov)*t_alpha-0.5*(cov-sigmasq)*(cov*(T)^(-2)*(y_dt_1'*y_dt_1))^(-1/2);
% end

function [p_zalpha,p_zt]=CV(q,dt,zalpha,zt,type,b)
load('critical_value_of_PP(fb).mat','cv_zas_bt','cv_zts_bt','cv_zas_pr','cv_zts_pr','cv_zas_qs','cv_zts_qs');
sigLevels=[0.001,(0.005:0.005:0.1),(0.125:0.025:0.8),(0.8005:0.005:0.995),0.999];
bs=(0.02:0.02:1);
switch type
    case 'bt'
        if q==-1
            CVTableZalpha=squeeze(cv_zas_bt(:,1,:));
            CVTableZt=squeeze(cv_zts_bt(:,1,:));
        elseif q==0 && dt==1
            CVTableZalpha=squeeze(cv_zas_bt(:,2,:));
            CVTableZt=squeeze(cv_zts_bt(:,2,:));
        elseif q==0 && dt==2
            CVTableZalpha=squeeze(cv_zas_bt(:,3,:));
            CVTableZt=squeeze(cv_zts_bt(:,3,:));
        elseif q==1 && dt==1
            CVTableZalpha=squeeze(cv_zas_bt(:,4,:));
            CVTableZt=squeeze(cv_zts_bt(:,4,:));
        elseif q==1 && dt==2
            CVTableZalpha=squeeze(cv_zas_bt(:,5,:));
            CVTableZt=squeeze(cv_zts_bt(:,5,:));
        end
    case 'pr'
        if q==-1
            CVTableZalpha=squeeze(cv_zas_pr(:,1,:));
            CVTableZt=squeeze(cv_zts_pr(:,1,:));
        elseif q==0 && dt==1
            CVTableZalpha=squeeze(cv_zas_pr(:,2,:));
            CVTableZt=squeeze(cv_zts_pr(:,2,:));
        elseif q==0 && dt==2
            CVTableZalpha=squeeze(cv_zas_pr(:,3,:));
            CVTableZt=squeeze(cv_zts_pr(:,3,:));
        elseif q==1 && dt==1
            CVTableZalpha=squeeze(cv_zas_pr(:,4,:));
            CVTableZt=squeeze(cv_zts_pr(:,4,:));
        elseif q==1 && dt==2
            CVTableZalpha=squeeze(cv_zas_pr(:,5,:));
            CVTableZt=squeeze(cv_zts_pr(:,5,:));
        end
    case 'qs'
        if q==-1
            CVTableZalpha=squeeze(cv_zas_qs(:,1,:));
            CVTableZt=squeeze(cv_zts_qs(:,1,:));
        elseif q==0 && dt==1
            CVTableZalpha=squeeze(cv_zas_qs(:,2,:));
            CVTableZt=squeeze(cv_zts_qs(:,2,:));
        elseif q==0 && dt==2
            CVTableZalpha=squeeze(cv_zas_qs(:,3,:));
            CVTableZt=squeeze(cv_zts_qs(:,3,:));
        elseif q==1 && dt==1
            CVTableZalpha=squeeze(cv_zas_qs(:,4,:));
            CVTableZt=squeeze(cv_zts_qs(:,4,:));
        elseif q==1 && dt==2
            CVTableZalpha=squeeze(cv_zas_qs(:,5,:));
            CVTableZt=squeeze(cv_zts_qs(:,5,:));
        end
end

CVTableZalphaRowb = interp2(sigLevels,bs,CVTableZalpha,sigLevels,b,'linear');
if zalpha <= CVTableZalphaRowb(1)
    
    p_zalpha = sigLevels(1);
    
elseif zalpha >= CVTableZalphaRowb(end)
    
    p_zalpha = sigLevels(end);
    
else
    
    p_zalpha = interp1(CVTableZalphaRowb,sigLevels,zalpha,'linear');
    
end

CVTableZtRowb = interp2(sigLevels,bs,CVTableZt,sigLevels,b,'linear');
if zt <= CVTableZtRowb(1)
    
    p_zt = sigLevels(1);
    
elseif zt >= CVTableZtRowb(end)
    
    p_zt = sigLevels(end);
    
else
    
    p_zt = interp1(CVTableZtRowb,sigLevels,zt,'linear');
    
end

end


