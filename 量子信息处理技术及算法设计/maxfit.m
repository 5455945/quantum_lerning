function [ max, c] = maxfit( fit )
% ����Ⱥ����Ӧ�����ĸ���
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

