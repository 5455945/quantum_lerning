function [ measurepop ] = measure( pop )
% measure(pop) ����Ⱥ�����������������������ƴ�
    global popsize;
    global lchrom;
    measurepop = zeros(popsize, lchrom);
    for i = 1:popsize
        for j = 1:lchrom
            x = rand(); %�������һ��0��1֮���һ����
            if x > pop(i, j, 1)^2
                measurepop(i, j) = 1;
            else
                measurepop(i, j) = 0;
            end
        end
    end
end

