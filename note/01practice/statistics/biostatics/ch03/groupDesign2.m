clear
x_A = [50.00,47.00,42.00,43.00,39.00,51.00,43.00,38.00,44.00,37.00];
x_B = [36.00,38.00,37.00,38.00,36.00,39.00,37.00,35.00,33.00,37.00];
x =[x_A,x_B];
group = [1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2];

% 方差齐性检验，即检验两组样本的总体方差是否相同
[p,stats] = vartestn(x',group','TestType','LeveneAbsolute'); %调用vartestn函数进行方差齐次性检验
%[p,stats] = vartestn(x',group','TestType','LeveneAbsolute','Display','off');

disp('Independent t-test:');
disp(['Levene’s test: p = ',num2str(p,'%0.2f')]);%方差检验方法：Levene检验

if p < 0.05
    disp('Equal variances not assumed') %方差不相同
    [h4,p4,ci4,stats4]=ttest2(x_A,x_B,'Vartype','unequal');
else
    disp('Equal variances assumed'); %方差相同
    [h4,p4,ci4,stats4]=ttest2(x_A,x_B);
end
disp(['t = ',num2str(stats4.tstat,'%0.2f')]);
disp(['df = ',num2str(stats4.df,'%0.2f')]);
disp(['p = ',num2str(p4,'%0.2f')]);