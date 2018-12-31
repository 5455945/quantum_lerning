# 《量子信息处理及算法设计》
  周日贵 科学出版社 2013.06

  这本书中除grover对比相关的代码，都收录了；
  除QSS中的t函数其它都可以在matlab r2016b上面跑通。
  
## 外积 op.m
op.m
```
function [outerproduct] = op(x, y)
    [m, n] = size(x);
	[p, q] = size(y);
    X = sprintf('%d %d %d %d', m, n, p, q);
    disp(X);
	if n ~= 1 || q ~= 1 || m ~= p
	    error('请正确输入做外积的2个向量')
    else
	    outerproduct = x*y';
    end
end
```

## 张量积 kron
```
a = [1 2 3];
b = [1 2 3; 4 5 6];
c = kron(a', b)
```

## Hadamard门 HS.m
```
function [Hadamards] = HS(x, n)  % n次的Hadamard门
    [a, b] = size(x);
	I = [1 1; 1 -1]/sqrt(2);
	if n > 1 && (n == fix(n))
	    temp = I;
	    for i = 1: n-1
		    temp = t(temp, I);
		end
		I = temp;
		[c, d] = size(I);
		if d == a
		    temp = I * x;
			Hadamards = temp;
	    else
		    error('Hadamard门不能作用此量子态');
		end
	else
	    error('请正确输入参数,作用次数必须大于等于2');
	end
end
```

## 受控非门 CNOT.m
```
function [controlnot] = CNOT(x)   % 多量子受控非门
    y = x;
	[a, b] = size(x);
	if b == 2
	    if (y(1) == '0' || y(1) == '1') && (y(2) == '0' || y(2) == '1')  % 输入判定
		    c = num2str(xor(str2num(y(1)), str2num(y(2))));
			controlnot = [y(1) c];
		else
		    error('请输入两位0、1字符组合');
		end
    else
	    error('受控非门不能作用此量子态');
	end
end
```

## 多量子比特量子态表示 QSS.m
```
% 这个函数里面没有找到t函数的定义
function quantumss = QSS(x) % 多量子比特
y = x;
[a, b] = size(x);
temp = 0;
for i = 1:b
    if (y(i) ~= '0') && (y(i) ~= '1')
        temp = 1;
    end
end

if temp == 0
    I = Q(str2num(y(1)));
    if b == 1
        quantumss = I;
    else
        for i = 1:b-1
            I = t(I, Q(str2num(y(i+1))));
        end
        quantumss = I;
    end
else
    SS = class(x);
    if SS(:, 1:4) =='doub'
        x = num2str(x);
        y = x;
        [a, b] = size(x);
        flag = 0;
        for i = 1:b
            if (y(i) ~= '0') && (y(i) ~= '1')
                flag = 1;
            end
        end
        if flag == 0
            I = Q(str2num(y(1)));
            if b == 1
                quantumss = I;
            else
                for i = 1:b-1
                    I = t(I, Q(str2num(y(i+1))));
                end
                quantumss = I;
            end
        else
            temp = str2num(x);
            temp = dec2bin(temp);
            [c, d] = size(temp);
            I = Q(str2num(temp(1)));
            for i = 1:d-1
                I = t(I, Q(str2num(temp(i+1))));
            end
            quantumss = I;
        end
    else
        temp = str2num(x);
        temp = dec2bin(temp);
        [c, d] = size(temp);
        I = Q(str2num(temp(1)));
        for i = 1:d-1
            I = t(I, Q(str2num(temp(i+1))));
        end
        quantumss = I;
    end
end
end
```

## 单量子比特量子态表示 Q.m
```
function [ quantum] = Q( x )
    i = x;
    if i == 0
        quantum = [1;0];
    else
        if i == 1
            quantum = [0;1];
        else
            disp('x 必须为0或1');
        end
    end
end
```

## n次张量积 ts.m
```
function [tensors] = ts(x, n)  % n次的张量积
    I = x;
	if n > 1 && (n == fix(n))
	    temp = I;
        for i = 1 : n - 1
            temp = t(temp, I);
        end
        I = temp;
        tensors = I;
    else
        disp('请输入整数次')
    end
end
% A = [1, 1, 2, 3; 4, 5, 6, 6];
% ts(A, 2)
```

## Bloch球面 bl.m
```
function [ bloch ] = bl( x ) % Bloch球面
    a = 2 * acos(x(1));
    if a == 0
        b = 0;
    else
        b = acos(real(x(2))/sin(a/2));
    end
    t = linspace(0, pi * 2, 37);
    p = linspace(0, pi, 19);
    [t, p] = meshgrid(t, p);
    r = 1;
    x = r.* sin(p).* cos(t);
    y = r.* sin(p).* sin(t);
    z = r.* cos(p);
    zz = z; zzz = z;
    % zz(p < pi/6) = nan; zzz(p > pi/6) = nan;
    surf(x, y, zz, 'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none');
    hold on
    % surf(x, y, zzz, 'facecolor', [1 1 0], 'edgecolor', 'none');
    camlight; lighting gouraud; alpha .5; axis equal;
    pp = linspace(0, 2 * pi, 19);
    h = linspace(0, 0, 19);
    plot3(sin(pp), cos(pp), h, '--')
    quiver3([0, 0], [0, 0], [0 0 ], [0 sin(a)*cos(b)], [0 sin(a)*sin(b)], [0 cos(a)])
end

% bl(QSS('0'))
% bl(QSS(1))
% c = sqrt(2)/2;
% bl(c*QSS(0) + c*QSS(1));
% bl(c*QSS(0) + exp(i*pi/6)*c*QSS(1));
```

## 判定是否具有酉性 u.m
```
function [ unitarity ] = u( x ) % 判定是否具有酉性
    [m, n] = size(x);
    if m ~= n
        disp('次矩阵不具有酉性');
    else
        I = x*x';
        J = eye(m, n);
        if det(J-I) < (10^-30) && det(J-I) > (-10^-30)
            disp('次矩阵具有酉性');
        else
            disp('次矩阵不具有酉性');
        end
    end
end

% I = [1 1; 1 -1]/sqrt(2);
% u(I)
```

## 判定是否是量子叠加态 SUP.m
```
function [ superposition ] = SUP( x ) % 判定是否是量子叠加态
    [a, b] = size(x);
    s = x(1).^2;
    for i = 1:a - 1
        s = s + x(i+1).^2;
    end
    if x(1) == 1 && x(a) == 0
        disp('次状态是量子态');
    elseif x(1) == 0 && x(a) == 1
        disp('次状态是量子激发态');
    elseif (1 - s) <= (10^-30) && (1-s) >= (-10^-30)
        disp('次状态是量子叠加态');
    else
        disp('次状态不是量子叠加态');
    end
end

% SUP(QSS('0'))
% SUP(QSS('1'))
% SUP(0.70/2*QSS('0') + 0.70*QSS('1'))
% SUP(sqrt(2)/2*QSS('0') + sqrt(2)/2*QSS('1'))
% SUP(0.5*QSS('00') + 0.5*QSS('01') + 0.5*QSS('10') + 0.5*QSS('11'))
```

## Toffoli门 TG.m
```
function [ Toffoli ] = TG( x ) % Toffoli门
    y = x;
    [a, b] = size(x);
    if b == 3
        if (y(1) == '0' || y(1) == '1') && (y(2) == '0' || y(2) == '1') && (y(3) == '0' || y(3) == '1')
            c = num2str(xor(str2num(y(3)), and (str2num(y(1)), str2num(y(2)))));
            Toffoli = [y(1) y(2) c];
        else
            disp('请输入三位0、1字符组合')
        end
    else
        disp('受控非门不能作用此量子态')
    end
end

% a = ['000'; '001'; '010'; '011'; '100'; '101'; '110'; '111'];
% for i = 1: 8
%     b = TG(a(i, :));
%     sprintf('%c %c %c 通过TG门后的值为 %c %c %c \n', a(i, 1), a(i, 2), a(i, 3),
%     b(1), b(2), b(3))
% end
```

## Fredkin门 FRG.m
```
function [ Fredkin ] = FRG( x ) % Fredkin 门
    y = x;
    [a, b] = size(x);
    if b == 3
        if (y(1) == '0' || y(1) == '1') && (y(2) == '0' || y(2) == '1') && (y(3) == '0' || y(3) == '1')
            if y(1) == '0'
                c = y(2);
                d = y(3);
            else
                c = y(3);
                d = y(2);
            end
            Fredkin = [y(1) c d];
        else
            disp('请输入三位0、1字符组合')
        end
    else
        disp('受控非门不能作用此量子态')
    end
end

% a = ['000'; '001'; '010'; '011'; '100'; '101'; '110'; '111'];
% for i = 1: 8
%     b = FRG(a(i, :));
%     sprintf('%c %c %c 通过FRG门后的值为 %c %c %c \n', a(i, 1), a(i, 2), a(i, 3), b(1), b(2), b(3))
% end
```

## 交换门 swap.m
```
function [ swap ] = swap( x, m, n)
    [a, b] = size(x);
    if (b < n) || (b < m) || (m == n) || n ~= fix(n) || m ~= fix(m)
        disp('请正确输入要交换的位')
    else
        temp = x(m);
        x(m) = x(n);
        x(n) = temp;
        swap = x;
    end
end

% swap('0111', 1, 2)
% swap('0111', 2, 1)
% swap('0111', 2, 1.2)
% swap('0111', 2, 5)
```


## 量子遗传算法
### 主运行脚本 QGA.m
```
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
```

### 种群初始化 intipop.m
```
function [ b ] = intipop(  )
% QGA算法产生初始化种群
    global popsize;
    global lchrom;
    b = zeros(popsize, lchrom, 2);
    for i = 1:popsize
        for j = 1:lchrom
            for k = 1:2
                b(i, j, k) = 1/sqrt(2);
            end
        end
    end
end
```

### cross.m
```
function [ newpop ] = cross( oldpop )
% cross
    global popsize;
    global lchrom;
    newpop = zeros(popsize, lchrom, 2);
    for j = 1:popsize
        for i = 1:popsize
            if i - j + 1 > 0
                c = i - j + 1;
            else
                c = i - j + 1 + popsize;
            end
            for k = 1:2
                newpop(i, j, k) = oldpop(c, j, k);
            end
        end
    end
    for j = (popsize + 1):lchrom
        for i = 1:popsize
            if i - j + 1 + popsize > 0
                c = i - j + 1 + popsize;
            else
                c = i - j + 1 + 2*popsize;
            end
            for k = 1:2
                newpop(i, j, k) = oldpop(c, j, k);
            end
        end
    end
end
```

### 测量函数 measure.m
```
function [ measurepop ] = measure( pop )
% measure(pop) 对种群进行塌缩测量，产生二进制串
    global popsize;
    global lchrom;
    measurepop = zeros(popsize, lchrom);
    for i = 1:popsize
        for j = 1:lchrom
            x = rand(); %随机产生一个0与1之间的一个数
            if x > pop(i, j, 1)^2
                measurepop(i, j) = 1;
            else
                measurepop(i, j) = 0;
            end
        end
    end
end
```

### 计算适应度 objfunc.m
```
function [ incode, fitness ] = objfunc( mepop )
% 对测量结果二进制数列进行解码,并计算适应度
% incode:种群的解码数列
% fitness:种群的适应度
    global popsize
    global lchrom;
    fitness = zeros(1, popsize);
    incode = zeros(1, popsize);
    % 解码
    for i = 1:popsize
        sum = 0;
        for j = 1:lchrom
            sum = sum + mepop(i,j)*2^(4-j);
        end
        incode(i) = sum;
    end
    %计算适应度
    for i = 1:popsize
        fitness(i) = exp(-0.001*incode(i))*(cos(0.8*incode(i)))^2;
    end
end
```

### 适应度最大个体 maxfit.m
```
function [ max, c] = maxfit( fit )
% 求种群中适应度最大的个体
    global popsize;
    global lchrom;
    max = 0;
    for i = 1:popsize
        if fit(i) > max
            max = fit(i);
            c = i;
        end
    end
end
```

### generation.m
```
function [ newpop ] = generation( oldpop, f, L )
% f为目标值, L为目标值对应的个体
    global popsize;
    global lchrom;
    %ang = angle(oldpop, f, L);
    ang = angle(oldpop);
    for i = 1:popsize
        for j = 1:lchrom
            H(:,:) = oldpop(i,j,:);
            a = ang(i,j);
            newpop(i,j,:) = [cos(a), -sin(a); sin(a), cos(a)] * H;
        end
    end
end
```


## Grover算法与pi/2相移算法的成功概率对比 GroverTest01.mlx
```
% 01 标准Grover算法与pi/2相移算法的成功概率对比
% 
% 标记态与状态总数比值
t = 0.001:0.01:1;
% 标准的Grover算法迭代步数
r = round(acos(sqrt(t))./(2*asin(sqrt(t))));
% 标准Grover算法成功概率
P = sin((2*r+1).*asin(sqrt(t))).^2;
% 当t > 1/3 且旋转相位位pi/2时，经过一步迭代后算法成功的概率
P1 = 4*t.^3 - 8*t.^2 + 5*t;
plot(t, P, t, P1);
xlabel('标记态与状态总数比值')
ylabel('成功概率')
legend('P', 'P1')
grid on
```
