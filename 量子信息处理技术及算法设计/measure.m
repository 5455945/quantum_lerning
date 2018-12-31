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

