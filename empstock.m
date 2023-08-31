clear
clc

rng(0)

data=readtable('stock index.xlsx');
time=data.time;
t=datenum(time);
HS300=data.HS300; % CSI 300 index

X=HS300;


plot(t,X);
xlim([t(1,1),t(end,1)]);
datetick('x','dd/mm/yyyy','keeplimits','keepticks');



dt=3;pp=4;
q=1;
ic='aic';
B=499;
b=0.02;
type='bt';
[zalpha,p_zalpha,zt,p_zt]=ztest(X,dt,q,pp,ic,B,b,type);
dtindex=detrend(X,1,3,0);
r=dtindex;
model=arima(1,1,1);
model.Constant=0;
c=estimate(model,r);