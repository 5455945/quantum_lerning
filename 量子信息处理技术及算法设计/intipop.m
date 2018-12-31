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

