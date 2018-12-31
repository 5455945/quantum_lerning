clear
tic;
clf
clc
global popsize;    % ��Ⱥ��С
global lchrom;     % Ⱦɫ�峤��
popsize = 10;
lchrom = 20;
oldpop = intipop;  %������ʼ��Ⱥ,���ѳ�ʼ��Ⱥ��ֵ�� oldpop
for i = 1:60
    oldpop = cross(oldpop)
    % ��Խ��棬�����ǰ�������������Ĳ��ּ����滻����������¸���Ĳ�����
    % �Ŵ��㷨����������õľ��ǽ�������
    mepop = measure(oldpop)
    [incode, fitness] = objfunc(mepop);
    % �Բ���������������н��н���,��������Ӧ�� incode ��Ⱥ�Ľ�������
    [max, c] = maxfit(fitness); %����Ⱥ����Ӧ�����ĸ���
    if (i == 1 | max > f)
        f = max;
        L = mepop(c, :);
        % mepop = measure(pop); measure ����Ⱥ�����������������������ƴ�
        x(i) = incode(c);
    else
        x(i) = x(i-1); % ���б���
    end
    maxnum(i) = f;
    newpop = generation(oldpop, f, L); % fΪĿ��ֵ��LΪĿ��ֵ��Ӧ�ĸ���
    old = newpop;
end
plot(maxnum, 'linewidth', 1.5),
title('�����Ŵ��㷨1�������Ž�'),
xlabel('��������'),
ylabel('����ֵ');
for i = 1:5:25
    maxnum(i)'
    x(i),
end
