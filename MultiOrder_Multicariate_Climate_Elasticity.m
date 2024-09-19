%% 需要更改的参数、路径
% 数据文件路径
path = 'YourOwnDataFilePath';
data = xlsread(path);

%% 主程序
numStation = data(1,1);

%% 定义变量
alpha = 0.05;

% 根据站点个数循环
for i = 1:numStation
    cf{i} = data(1:end,i*8);
    snow{i} = data(1:end,i*8-1);
    Pr{i} = data(1:end,i*8-2);
    T{i} = data(1:end,i*8-3);
    QS{i} = data(1:end,i*8-4);
    Q{i} = data(1:end,i*8-5);
    year{i} = data(1:end,i*8-6);

    % 所有变量整合 删除NaN年份
    allData_deleteNan{i} = [year{i},Q{i},QS{i},T{i},Pr{i},snow{i},cf{i}];
    allData_deleteNan{i} = allData_deleteNan{i}(all(~isnan(allData_deleteNan{i}),2),:);

    cf{i} = allData_deleteNan{i}(:,7);
    snow{i} = allData_deleteNan{i}(:,6);
    Pr{i} = allData_deleteNan{i}(:,5);
    T{i} = allData_deleteNan{i}(:,4);
    QS{i} = allData_deleteNan{i}(:,3);
    Q{i} = allData_deleteNan{i}(:,2);
    year{i} = allData_deleteNan{i}(:,1);

    num_year = size(year{i},1);

    % 三年滑动平均
    for j = 1:num_year-2
        cf_3year_average{i}(j,1) = nanmean(cf{i}(j:j+2),"all");
        snow_3year_average{i}(j,1) = nanmean(snow{i}(j:j+2),"all");
        Pr_3year_average{i}(j,1) = nanmean(Pr{i}(j:j+2),"all");
        T_3year_average{i}(j,1) = nanmean(T{i}(j:j+2),"all");
        QS_3year_average{i}(j,1) = nanmean(QS{i}(j:j+2),"all");
        Q_3year_average{i}(j,1) = nanmean(Q{i}(j:j+2),"all");
    end

    % 各变量所有年平均
    cf_average{i} = nanmean(cf_3year_average{i},"all");
    snow_average{i} = nanmean(snow_3year_average{i},"all");
    Pr_average{i} = nanmean(Pr_3year_average{i},"all");
    T_average{i} = nanmean(T_3year_average{i},"all");
    QS_average{i} = nanmean(QS_3year_average{i},"all");
    Q_average{i} = nanmean(Q_3year_average{i},"all");

    % 各变量所有年最大 最小
    cf_max{i} = max(cf_3year_average{i});
    cf_min{i} = min(cf_3year_average{i});
    snow_max{i} = max(snow_3year_average{i});
    snow_min{i} = min(snow_3year_average{i});
    Pr_max{i} = max(Pr_3year_average{i});
    Pr_min{i} = min(Pr_3year_average{i});
    T_max{i} = max(T_3year_average{i});
    T_min{i} = min(T_3year_average{i});
    QS_max{i} = max(QS_3year_average{i});
    QS_min{i} = min(QS_3year_average{i});
    Q_max{i} = max(Q_3year_average{i});
    Q_min{i} = min(Q_3year_average{i});


    % 回归变量 = (X-Xmean)/(Xmax-Xmin)
    cf_3year_normalization_change{i} = (cf_3year_average{i}-cf_average{i}) / (cf_max{i}-cf_min{i});
    snow_3year_normalization_change{i} = (snow_3year_average{i}-snow_average{i})/(snow_max{i}-snow_min{i});
    Pr_3year_normalization_change{i} = (Pr_3year_average{i}-Pr_average{i})/(Pr_max{i}-Pr_min{i});
    T_3year_normalization_change{i} = (T_3year_average{i}-T_average{i})/(T_max{i}-T_min{i});
    QS_3year_normalization_change{i} = (QS_3year_average{i}-QS_average{i})/(QS_max{i}-QS_min{i});
    Q_3year_normalization_change{i} = (Q_3year_average{i}-Q_average{i})/(Q_max{i}-Q_min{i});

    %% 所有路径回归计算
    %% 1 Qs_P_Q_T_V
    y = QS_3year_normalization_change{i};
    x1 = Pr_3year_normalization_change{i};
    x2 = Q_3year_normalization_change{i};
    x3 = T_3year_normalization_change{i};
    x4 = cf_3year_normalization_change{i};

    % 去除NaN
    delete_nan =[y,x1,x2,x3,x4];
    delete_nan = delete_nan(all(~isnan(delete_nan),2),:);

    % 回归变量
    % 创建表格，定义名称
    regress_x = table(delete_nan);

    regress_y = table(delete_nan(:,1),'VariableNames',{'QS'});
    regress_x = table(delete_nan(:,2),delete_nan(:,3),delete_nan(:,4),delete_nan(:,5), ...
        delete_nan(:,2).^2,delete_nan(:,3).^2,delete_nan(:,4).^2,delete_nan(:,5).^2, ...
        delete_nan(:,2).^3,delete_nan(:,3).^3,delete_nan(:,4).^3,delete_nan(:,5).^3, ...
        delete_nan(:,2).^4,delete_nan(:,3).^4,delete_nan(:,4).^4,delete_nan(:,5).^4, ...
        delete_nan(:,2).^5,delete_nan(:,3).^5,delete_nan(:,4).^5,delete_nan(:,5).^5, ...
        'VariableNames',{'P','Q','T','V', ...
        'P_2','Q_2','T_2','V_2', ...
        'P_3','Q_3','T_3','V_3', ...
        'P_4','Q_4','T_4','V_4', ...
        'P_5','Q_5','T_5','V_5'});

    % 回归计算
    mdl_Qs_P_Q_T_V_stepPathway{i} = stepwiseFoward(regress_x,regress_y,alpha);

    %% 2 Q_P_T_V
    y = Q_3year_normalization_change{i};
    x1 = Pr_3year_normalization_change{i};
    x2 = T_3year_normalization_change{i};
    x3 = cf_3year_normalization_change{i};

    % 去除NaN
    delete_nan =[y,x1,x2,x3];
    delete_nan = delete_nan(all(~isnan(delete_nan),2),:);

    % 回归变量
    % 创建表格，定义名称
    regress_x = table(delete_nan);

    regress_y = table(delete_nan(:,1),'VariableNames',{'Q'});
    regress_x = table(delete_nan(:,2),delete_nan(:,3),delete_nan(:,4), ...
        delete_nan(:,2).^2,delete_nan(:,3).^2,delete_nan(:,4).^2, ...
        delete_nan(:,2).^3,delete_nan(:,3).^3,delete_nan(:,4).^3, ...
        delete_nan(:,2).^4,delete_nan(:,3).^4,delete_nan(:,4).^4, ...
        delete_nan(:,2).^5,delete_nan(:,3).^5,delete_nan(:,4).^5, ...
        'VariableNames',{'P','T','V', ...
        'P_2','T_2','V_2', ...
        'P_3','T_3','V_3', ...
        'P_4','T_4','V_4', ...
        'P_5','T_5','V_5'});

    % 回归计算
    mdl_Q_P_T_V_stepPathway{i} = stepwiseFoward(regress_x,regress_y,alpha);

    %% 3 snow_T
    y = snow_3year_normalization_change{i};
    x1 = T_3year_normalization_change{i};

    % 去除NaN
    delete_nan =[y,x1];
    delete_nan = delete_nan(all(~isnan(delete_nan),2),:);

    % 回归变量
    % 创建表格，定义名称
    regress_x = table(delete_nan);

    regress_y = table(delete_nan(:,1),'VariableNames',{'SNOW'});
    regress_x = table(delete_nan(:,2),delete_nan(:,2).^2,delete_nan(:,2).^3,delete_nan(:,2).^4,delete_nan(:,2).^5, ...
        'VariableNames',{'T','T_2','T_3','T_4','T_5'});

    % 回归计算
    mdl_SNOW_T_stepPathway{i} = stepwiseFoward(regress_x,regress_y,alpha);
  
    %% 结果整合
    regress_Result_stepPathway{i}.Qs_P_Q_T_V = mdl_Qs_P_Q_T_V_stepPathway{i};
    regress_Result_stepPathway{i}.Q_P_T_V = mdl_Q_P_T_V_stepPathway{i};
    regress_Result_stepPathway{i}.SNOW_T = mdl_SNOW_T_stepPathway{i};
   
end

%% 逐步回归 "前进法" 函数
function [mdl] = stepwiseFoward(regress_x,regress_y,alpha)
% 逐步回归 "前进法"
regress_x_initial_changing = regress_x;
regress_x_final = table();

% 相当于第一个循环
num_regress_x_initial_changing = size(regress_x_initial_changing,2);
num_regress_x_final = size(regress_x_final,2);

% 每个变量单独与Y进行回归，找出t检验p值最小的变量，添加
for j = 1:num_regress_x_initial_changing
    regress_x_initial_name_single{j} = regress_x_initial_changing.Properties.VariableNames{j};
    regress_x_middle = [regress_x_final,regress_x_initial_changing(:,regress_x_initial_name_single{j})];

    mdl_middle{j} = fitlm([regress_x_middle,regress_y]);
    pValue_regress_add_x_middle(j) = mdl_middle{j}.Coefficients{num_regress_x_final+2,4};
    pValue_regress_all_x_middle{j} = mdl_middle{j}.Coefficients{2:end,4};
end

% 确定要增加的变量名
order_pValueMin_regress_x_initial_changing = find(pValue_regress_add_x_middle == min(pValue_regress_add_x_middle));
regress_x_add_name = regress_x_initial_changing.Properties.VariableNames(order_pValueMin_regress_x_initial_changing);

% 更改变量
regress_x_initial_changing = removevars(regress_x_initial_changing,regress_x_add_name);
num_regress_x_initial_changing = size(regress_x_initial_changing,2);

% 正式循环
% 如果还有变量p值小于0.05
while min(pValue_regress_add_x_middle) < alpha && max(pValue_regress_all_x_middle{order_pValueMin_regress_x_initial_changing}) < alpha

    % 更新各变量
    regress_x_final = [regress_x_final regress_x(:,regress_x_add_name)];
    num_regress_x_final = size(regress_x_final,2);

    pValue_regress_add_x_middle = [];
    pValue_regress_all_x_middle = {};

    % 再次回归
    for j = 1:num_regress_x_initial_changing
        regress_x_initial_name_single{j} = regress_x_initial_changing.Properties.VariableNames{j};
        regress_x_middle = [regress_x_final,regress_x_initial_changing(:,regress_x_initial_name_single{j})];

        mdl_middle{j} = fitlm([regress_x_middle,regress_y]);
        pValue_regress_add_x_middle(j) = mdl_middle{j}.Coefficients{num_regress_x_final+2,4};
        pValue_regress_all_x_middle{j} = mdl_middle{j}.Coefficients{2:end,4};
    end

    % 确定要增加的变量名
    order_pValueMin_regress_x_initial_changing = find(pValue_regress_add_x_middle == min(pValue_regress_add_x_middle));
    regress_x_add_name = regress_x_initial_changing.Properties.VariableNames(order_pValueMin_regress_x_initial_changing);

    % 更改变量
    regress_x_initial_changing = removevars(regress_x_initial_changing,regress_x_add_name);
    num_regress_x_initial_changing = size(regress_x_initial_changing,2);

end

% 最终回归变量确定，进行回归
mdl = fitlm([regress_x_final,regress_y]);
end