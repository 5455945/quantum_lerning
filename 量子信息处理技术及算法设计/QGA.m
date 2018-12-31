clear
tic;
clf
clc
global popsize;    % 种群大小
global lchrom;     % 染色体长度
popsize = 10;
lchrom = 20;
oldpop = intipop;  %产生初始种群,并把初始种群赋值给 oldpop
for i = 1:60
    oldpop = cross(oldpop)
    % 配对交叉，交叉是把两个父代个体的部分加以替换重组而生成新个体的操作。
    % 遗传算法中起核心作用的就是交叉算子
    mepop = measure(oldpop)
    [incode, fitness] = objfunc(mepop);
    % 对测量结果二进制数列进行解码,并计算适应度 incode 种群的解码数列
    [max, c] = maxfit(fitness); %求种群中适应度最大的个体
    if (i == 1 | max > f)
        f = max;
        L = mepop(c, :);
        % mepop = measure(pop); measure 对种群进行塌缩测量，产生二进制串
        x(i) = incode(c);
    else
        x(i) = x(i-1); % 进行变体
    end
    maxnum(i) = f;
    newpop = generation(oldpop, f, L); % f为目标值，L为目标值对应的个体
    old = newpop;
end
plot(maxnum, 'linewidth', 1.5),
title('量子遗传算法1：求最优解'),
xlabel('进化代数'),
ylabel('最优值');
for i = 1:5:25
    maxnum(i)'
    x(i),
end
