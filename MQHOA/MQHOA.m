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
% -------������������ʼ
repeat = 10; % �����ظ��������
DIM = 300; % objective function's dimentions
ACCURACY = 0.000001; % ������㾫�ȣ�sigma_min
optimalNum = 15; % kֵ�Ĵ�С
minDomain = -10; maxDomain = 10; % objective function's domain
% -------���������������
% -------����ظ������������
gmeanV = zeros(1, repeat);
gbestV = zeros(1, repeat);
gfangcha = zeros(1, repeat);
% -------����ظ�����
for kk = 1: repeat
    % --������ʼ������ʼ
    funcV = zeros(1, optimalNum); % ���屣��k��������ĺ���ֵ����ʼ��Ϊ0
    samplePos = zeros(DIM, optimalNum); % ���������
    sigma = maxDomain - minDomain; % ��ʼ�߶�ΪĿ�꺯��������Ĵ�С
    stdPre = zeros(1, DIM); % �洢ÿһ��ά���ϵı�׼��
    stdNow = zeros(1, DIM); 
    % -- ����k������λ�õ�Ŀ�꺯��ֵ
    % ����k��DIMά����������꣬����ʼ��
    optimalSolution = unifrnd(minDomain, maxDomain, DIM, optimalNum);
    stdPre = std(optimalSolution, 1, 2); % �����ʼ��k��������ı�׼�������
    w = 0; % function evolution times
    for k = 1: optimalNum % �����Ž⺯��ֵ
        funcV(k) = func(optimalSolution(:,k), DIM);
        w = w + 1;
    end
    while 1 % M iteration begin
        while 1 % QHO iteration begin
            while 1 %�ȶ�������������ʼ
                for i = 1: optimalNum
                    % ����Box-Muller���������µ���̬�ֲ����������
                    theat = 2*pi*rand(DIM, 1);
                    R = sqrt(-2.0*log(rand(DIM, 1)));
                    gaosiRand = R.*cos(theat);
                    samplePos(:, i) = optimalSolution(:, i) + sigma*gaosiRand;
                    % ����Box-Muller���������µ���̬�ֲ�������
                    % ����Խ�������
                    for i1 = 1: DIM
                        if samplePos(i1, i) > maxDomain
                            samplePos(i1, i) = maxDomain;
                        end
                        if samplePos(i1, i) < minDomain
                            samplePos(i1, i) = minDomain;
                        end
                    end
                    sampleValue = func(samplePos(:, i), DIM); % ���i��������ĺ���ֵ
                    w = w + 1;
                    if sampleValue < funcV(i) % ���������ֵС�ڵ�ǰ�㺯��ֵ�����滻
                        funcV(i) = sampleValue;
                        optimalSolution(:, i) = samplePos(:, i);
                    end
                end
                stdNow = std(optimalSolution, 1, 2); % �½��׼��
                c = stdPre - stdNow; % �������Ž��׼��Ĳ�ֵ
                stdPre = stdNow;
                if max(abs(c)) <= (sigma) % �������ά�еı�׼������ֵ��С��Sigma��������ڸ��ܼ��������ȶ�
                    break;
                end
            end % �ȶ���������������
            % -- ���´��뽫k���������о�ֵ��λ�úϲ������ֵ��λ���У�ʵ���ܼ����½�
            meanV = mean(funcV);
            meanPos = mean(optimalSolution, 2); % ȡ��ƽ������
            [v_max, index_max] = max(funcV); % ȡ�����ֵ�����index_max
            optimalSolution(:, index_max) = meanPos; % ��ƽ�������滻���ֵ��Ӧ����
            funcV(index_max) = func(meanPos, DIM);
            w = w + 1;
            stdPre = std(optimalSolution, 1, 2); % �½��׼��
            % -- �ܼ��½�����
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
         