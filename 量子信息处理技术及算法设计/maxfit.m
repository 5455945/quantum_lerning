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

