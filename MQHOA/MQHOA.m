% ---------------------------------------------------
% MQHOA algorithm for sigle model objective function
% Matlab Code by Wang Peng (Dec.21, 2015)
% This program must include func.m
% This is global minimum program.
% All rights reserved by Parallel Computing Lab.
% ---------------------------------------------------
clear all;
feature jit off; % shut down the accelerator
format long;
global DIM;
% -------参数定义区域开始
repeat = 10; % 定义重复计算次数
DIM = 300; % objective function's dimentions
ACCURACY = 0.000001; % 定义计算精度，sigma_min
optimalNum = 15; % k值的大小
minDomain = -10; maxDomain = 10; % objective function's domain
% -------参数定义区域结束
% -------多次重复计算所需参数
gmeanV = zeros(1, repeat);
gbestV = zeros(1, repeat);
gfangcha = zeros(1, repeat);
% -------多次重复计算
for kk = 1: repeat
    % --参数初始化区域开始
    funcV = zeros(1, optimalNum); % 定义保存k个采样点的函数值，初始化为0
    samplePos = zeros(DIM, optimalNum); % 采样点矩阵
    sigma = maxDomain - minDomain; % 初始尺度为目标函数定义域的大小
    stdPre = zeros(1, DIM); % 存储每一个维度上的标准差
    stdNow = zeros(1, DIM); 
    % -- 计算k个采样位置的目标函数值
    % 定义k个DIM维采样点的坐标，并初始化
    optimalSolution = unifrnd(minDomain, maxDomain, DIM, optimalNum);
    stdPre = std(optimalSolution, 1, 2); % 计算初始化k个采样点的标准差，按行求
    w = 0; % function evolution times
    for k = 1: optimalNum % 求最优解函数值
        funcV(k) = func(optimalSolution(:,k), DIM);
        w = w + 1;
    end
    while 1 % M iteration begin
        while 1 % QHO iteration begin
            while 1 %稳定性收敛迭代开始
                for i = 1: optimalNum
                    % 采用Box-Muller方法生成新的正态分布随机采样点
                    theat = 2*pi*rand(DIM, 1);
                    R = sqrt(-2.0*log(rand(DIM, 1)));
                    gaosiRand = R.*cos(theat);
                    samplePos(:, i) = optimalSolution(:, i) + sigma*gaosiRand;
                    % 采用Box-Muller方法生成新的正态分布采样点
                    % 处理越界采样点
                    for i1 = 1: DIM
                        if samplePos(i1, i) > maxDomain
                            samplePos(i1, i) = maxDomain;
                        end
                        if samplePos(i1, i) < minDomain
                            samplePos(i1, i) = minDomain;
                        end
                    end
                    sampleValue = func(samplePos(:, i), DIM); % 求第i个采样点的函数值
                    w = w + 1;
                    if sampleValue < funcV(i) % 如果采样点值小于当前点函数值，则替换
                        funcV(i) = sampleValue;
                        optimalSolution(:, i) = samplePos(:, i);
                    end
                end
                stdNow = std(optimalSolution, 1, 2); % 新解标准差
                c = stdPre - stdNow; % 两次最优解标准差的差值
                stdPre = stdNow;
                if max(abs(c)) <= (sigma) % 如果所有维中的标准差最大的值都小于Sigma，则表明在该能级采样已稳定
                    break;
                end
            end % 稳定性收敛迭代结束
            % -- 以下代码将k个采样点中均值的位置合并到最大值的位置中，实现能级的下降
            meanV = mean(funcV);
            meanPos = mean(optimalSolution, 2); % 取得平均坐标
            [v_max, index_max] = max(funcV); % 取得最大值的序号index_max
            optimalSolution(:, index_max) = meanPos; % 用平均坐标替换最大值对应坐标
            funcV(index_max) = func(meanPos, DIM);
            w = w + 1;
            stdPre = std(optimalSolution, 1, 2); % 新解标准差
            % -- 能级下降结束
            if max(stdPre) < sigma
                break;
            end
        end % QHO iteration end
        sigma = sigma/2.0;
        if sigma <= ACCURACY
            break;
        end
    end % M iteration end
    [global_min, index] = min(funcV);
    gmeanV(kk) = mean(funcV);
    gbestV(kk) = min(funcV);
    fprintf('No. %d run, DIM=%d, Global minimum function value is %e. function evolution = %d\n', kk, DIM, global_min, w);
end % MQHOA end
meanV = mean(gmeanV);
bestV = min(gbestV);
fangcha = std(gbestV);
fprintf('Repeat=%d, MeanValue=%d, BestValue=%d, Std=%d\n', repeat, meanV, bestV, fangcha);
         