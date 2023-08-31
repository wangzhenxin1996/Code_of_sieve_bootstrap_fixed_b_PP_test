
clear
clc

load('PowerCurveforSampleSizewithSizeAdjusted.mat','size_zas','size_zts')


Ts=(25:25:500);
%% plot for zalpha
folder_path = 'PowerCurvewithSampleSize\zalpha';

if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

for j=1:5
    for k=1:3
        for i=1:4
            data=squeeze(size_zas(j,(k-1)*4+i,:));
            plot(Ts,data)
            hold on
        end
        ylabel('Rejection rates');
        xlabel('T')
        legend('PP', 'PP(fb)1step','PP(fb)2step','PP^b(fb)', 'Location', 'southoutside','NumColumns',4);
        hold off
        img_name = ['PowerCurveforZALPHAdt_', num2str(k),'_ar_',num2str(j),'_withT'];
        print(fullfile(folder_path,img_name),'-depsc','-painters')
    end
end

%% plot for zalpha
folder_path = 'PowerCurvewithSampleSize\zt';


if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

for j=1:5
    for k=1:3
        for i=1:4
            data=squeeze(size_zts(j,(k-1)*4+i,:));
            plot(Ts,data)
            hold on
        end
        ylabel('Rejection rates');
        xlabel('T')
        legend('PP', 'PP(fb)1step','PP(fb)2step','PP^b(fb)', 'Location', 'southoutside','NumColumns',4);
        hold off
        img_name = ['PowerCurveforZTdt_', num2str(k),'_ar_',num2str(j),'_withT'];
        print(fullfile(folder_path,img_name),'-depsc','-painters')
    end
end
