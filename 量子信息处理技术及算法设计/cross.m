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

