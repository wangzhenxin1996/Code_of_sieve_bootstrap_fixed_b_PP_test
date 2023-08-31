clear;
clc;
rng(0)

% [re] represents the number of repetitions for MC simulation
% [T] represents the length of simulated samples
% [dt] represents the detrending method: '1' for one-step detrending;
%     '2' for two-step detrending, '3' for GLS detrending
% [q] represents the order of detrending: q=-1 for no-detrending; q=0 for demeaning;
%     q=1 for linear trend detrending
% [c] represents the degree of deviation from a unit root: a=1+c/T;
% [var] represents the variance of the disturbance term in constructing MC samples;
% [B] represents the number of bootstrap samples constructed using the bootstrap method;
% [delta] represents the coefficient vector of the trend term when constructing MC samples;
% [pp] represents the method for PP-test: "1" for traditional PP test,
%     which uses Newey-West for long-run variance estimation;
%     "2" for fixed-b PP unit root test, PP(fb), proposed by Vogelsang and Wagner (2013),
%     which replaces Newey-West long-run variance estimation with fixed-b HAR long-run variance estimation;
%     "4" for obtaining critical values of PP(fb) test using Sieve Bootstrap method, PP^b(fb);
% [b] represents the truncation parameter for fixed smoothing parameter HAR long-run variance estimation:
%     when "type" is "bt", "pr", or "qs", 'b' represents the truncation parameter;
%     when "type" is "os", it represents the length of the orthogonal sequence used;
% [type] represents the type of kernel or non-parametric method used when estimating
%     long-run variance-covariance matrix: "bt" for Bartlett kernel;
%     "pr" for Parzen kernel; "qs" for QuadraticSpectral kernel;
%     "os" for Orthonormal Series method.
%
% **Reference:**
% Vogelsang, T. J. and Wagner, M. (2013). A Fixed-b Perspective on the Phillips-Perron Unit Root Tests.
% Econometric Theory, 29(3):609-628.


re=2000;

cs=(-1.25:-1.25:-25);

var=1;
B=499;
delta=[0;0];
ic='aic';

b=0.02;
Ts=(25:25:500);
num_Ts=length(Ts);
qs=[-1 0 1];
rhos=[0, 0.5, -0.3, 0, 0;
    0,  0,  -0.2, 0, 0];
phis=[0, 0, 0, -0.5, 0.3;
    0, 0, 0,   0,  0.2];
pps=[1 2 2 4];
dts=[1 1 2 3];
type='bt';
size_zas=zeros(5,12,num_Ts);
size_zts=zeros(5,12,num_Ts);

load('CriticalValueforSizeAdjustedofPPandPP(fb).mat','cv_zas','cv_zts');

for h=1:num_Ts
    T=Ts(h);
    c=cs(h);
    for e=1:3 % Demeaned or Detrend
        q=qs(e);
        for f=1:5  % Constant or AR(1) or MA(1)
            rho=rhos(:,f)';
            phi=phis(:,f)';
            for g=1:4 % PP with OLS_detrend or PP(fb) with OLS detrend or PP^b(fb) with GLS detrend
                pp=pps(g);
                dt=dts(g);

                zalpha=zeros(re,1);
                p_zalpha=zeros(re,1);
                zt=zeros(re,1);
                p_zt=zeros(re,1);


                y=DGP(T,re,c,delta,rho,phi,var);
                parfor j=1:re
                    s = RandStream('mt19937ar','Seed',j);
                    RandStream.setGlobalStream(s);

                    [zalpha(j,1),p_zalpha(j,1),zt(j,1),p_zt(j,1)]=ztest(y(:,j),dt,q,pp,ic,B,b,type);
                end
                if pp==4
                    size_za=sum(p_zalpha(:,1)<=0.05)/re;
                    size_zt=sum(p_zt(:,1)<=0.05)/re;
                elseif pp==1 || pp==2
                    size_za=sum(zalpha(:,1)<=cv_zas(f,(e-1)*4+g,(h-1)*5+1))/re;
                    size_zt=sum(zt(:,1)<=cv_zts(f,(e-1)*4+g,(h-1)*5+1))/re;
                end
                size_zas(f,(e-1)*4+g,h)=size_za;
                size_zts(f,(e-1)*4+g,h)=size_zt;
            end
        end
    end
end

save('PowerCurveforSampleSizewithSizeAdjusted.mat','size_zas','size_zts');
