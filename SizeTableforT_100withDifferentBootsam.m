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

c=0;
delta=[0;0];
var=1;
B=499;
ic='aic';

b=0.02;

qs=[-1 0 1];

T=100;

Bs=[99 199 299 399 499];
rhos=[0, 0.5, -0.3, 0, 0;
    0,  0,  -0.2, 0, 0];
phis=[0, 0, 0, -0.5, 0.3;
    0, 0, 0,   0,  0.2];
pp=4;
dt=3;
type='bt';
size_zas=zeros(5,12);
size_zts=zeros(5,12);

% Nodetrended or Demeaned or Detrend
for e=1:3
    q=qs(e);
    % Constant or AR(1) or AR(2) or MA(1) or MA(2)
    for f=1:5  
        rho=rhos(:,f)';
        phi=phis(:,f)';
        % Different bootstrap number B  
        for g=1:5 
            B=Bs(g);

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

            size_za=sum(p_zalpha(:,1)<=0.05)/re;
            size_zt=sum(p_zt(:,1)<=0.05)/re;
            size_zas(f,(e-1)*5+g)=size_za;
            size_zts(f,(e-1)*5+g)=size_zt;
        end
    end
end
save('SizeTableforT_100withDifferentBootsam','size_zas','size_zts');
