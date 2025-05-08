clear
data = [506, 500, 495, 488, 504, 486, 505, 513, 521, 520, 512, 485];
% 计算平均值和标准差
mean_value = mean(data);
std_value = std(data);
% 计算t分布的临界值(t分布的反函数)
t_critical = tinv(0.975, length(data)-1);
% 计算置信区间
confidence_interval = t_critical * std_value / sqrt(length(data));
% 绘制数据点
plot(data, 'o');
hold on;
% 绘制置信区间
x = 1:length(data);
y = repmat(mean_value, 1, length(data));
upper_bound = y + confidence_interval
lower_bound = y - confidence_interval
fill([x, fliplr(x)], [upper_bound, fliplr(lower_bound)], 'b',
'FaceAlpha', 0.3);
% 添加图例
legend('数据点'
, '95%置信区间');
% 添加标题和轴标签
title('12个数据点的95%置信区间');
xlabel('数据点编号');
ylabel('数据值');

% 分别求出置信度为0.90，0.92，0.95，0.98的点估计和区间估计值
Pv=[0.9,0.92,0.95,0.98]; A=[];
disp(' 置信度 均值 标准差 均值下界 均值上界 标准差下界 标准差上界')
for i=1:length(Pv)
[muHat,sigmaHat,muCI,sigmaCI]=normfit(data,Pv(i));
A=[A; Pv(i),muHat,sigmaHat,muCI(1),muCI(2),sigmaCI(1),sigmaCI(2)];
end
disp(A)