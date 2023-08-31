
clear
clc

load('sizecurveforT_100withdifferentb.mat','size_zas','size_zts')


bs=(0.02:0.02:1);
y = 0.05 * ones(size(bs));

%% plot for zalpha
folder_path = 'SizeCurveT_100\zalpha';

if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

for j=1:5
    for k=1:3
        for i=1:4
            data=squeeze(size_zas(j,(k-1)*4+i,:));
            plot(bs,data)
            hold on
        end
        plot(bs, y, 'r--');
        ylabel('Rejection rates');
        xlabel('b')
        legend('PP', 'PP(fb)1step','PP(fb)2step','PP^b(fb)', 'Location', 'southoutside','NumColumns',4);
        hold off
        img_name = ['SizeCurveforZALPHAdt_', num2str(k),'_ar_',num2str(j),'_T_100'];
        print(fullfile(folder_path,img_name),'-depsc','-painters')
    end
end

%% plot for zalpha
folder_path = 'SizeCurveT_100\zt';

if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

for j=1:5
    for k=1:3
        for i=1:4
            data=squeeze(size_zts(j,(k-1)*4+i,:));
            plot(bs,data)
            hold on
        end
        plot(bs, y, 'r--');
        ylabel('Rejection rates');
        xlabel('b')
        legend('PP', 'PP(fb)1step','PP(fb)2step','PP^b(fb)', 'Location', 'southoutside','NumColumns',4);
        hold off
        img_name = ['SizeCurveforZTdt_', num2str(k),'_ar_',num2str(j),'_T_100'];
        print(fullfile(folder_path,img_name),'-depsc','-painters')
    end
end

