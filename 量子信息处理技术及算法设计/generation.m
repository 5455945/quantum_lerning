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

