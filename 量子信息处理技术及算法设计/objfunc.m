function [ incode, fitness ] = objfunc( mepop )
% �Բ���������������н��н���,��������Ӧ��
% incode:��Ⱥ�Ľ�������
% fitness:��Ⱥ����Ӧ��
    global popsize
    global lchrom;
    fitness = zeros(1, popsize);
    incode = zeros(1, popsize);
    % ����
    for i = 1:popsize
        sum = 0;
        for j = 1:lchrom
            sum = sum + mepop(i,j)*2^(4-j);
        end
        incode(i) = sum;
    end
    %������Ӧ��
    for i = 1:popsize
        fitness(i) = exp(-0.001*incode(i))*(cos(0.8*incode(i)))^2;
    end
end

