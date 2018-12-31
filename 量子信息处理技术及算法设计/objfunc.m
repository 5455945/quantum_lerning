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

