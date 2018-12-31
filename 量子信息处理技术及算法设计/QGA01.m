clear all;                                                   
close all;                                                   
%------------------------��������---------------------             
popsize = 100;        %��Ⱥ��ģ                              
vartotal = 2;         %����������һ��Ⱦɫ�������λ��          
shiftstep = 0.01*pi; %ת�ǲ���
Pm = ones(1,popsize)*0.05;%���ñ������
maxgen = 500;  %���õ�������
%------���鲿��--��ռ���Ż�������ȡֵ��Χ------------
var_range(1,1) = -100;    
var_range(1,2) = 100;    
var_range(2,1) = -100;  
var_range(2,2) = 100;     
%--------------------Ⱦɫ���ʼ��----------------------
%    ��ʼ����2*popsize��Ⱦɫ�壬����ÿ��Ⱦɫ������������
for i=1:1:vartotal
    
    fai(:,i)=2*pi*rand(popsize,1); 
    chrom(:,1,i)=cos(fai(:,i));
    chrom(:,2,i)=sin(fai(:,i));
   
    oldfai(:,i)=2*pi*rand(popsize,1); 
    oldchrom(:,1,i)=cos(oldfai(:,i));
    oldchrom(:,2,i)=sin(oldfai(:,i));
end
%--------------------��ռ�任-------------------------
for i=1:1:2 
    for j=1:1:vartotal
        chromx(:,i,j)=0.5*(var_range(j,2)*(1+chrom(:,i,j))+var_range(j,1)*(1-chrom(:,i,j)));
        oldchromx(:,i,j)=0.5*(var_range(j,2)*(1+oldchrom(:,i,j))+var_range(j,1)*(1-oldchrom(:,i,j)));
    end
end
%-----------������Ӧ��---��Ӧ�Ⱥ�����Shaffer's F6����-------------
for i=1:1:popsize
    for j=1:1:2
        x1=chromx(i,j,1);
        x2=chromx(i,j,2);
        fitness(i,j)=0.5-((sin(sqrt(x1^2+x2^2)))^2-0.5)/(1+0.001*(x1^2+x2^2))^2;
        x1=oldchromx(i,j,1);
        x2=oldchromx(i,j,2);
        oldfitness(i,j)=0.5-((sin(sqrt(x1^2+x2^2)))^2-0.5)/(1+0.001*(x1^2+x2^2))^2;
    end
end
%----------------------������Ž⼰��Ӧ�Ա���----------------------
[Bestf,Indexf]=sort(fitness,2);
[BestF,IndexF]=sort(Bestf,1);
gBestfit=BestF(popsize,2); 
gBestpop=IndexF(popsize,2); 
gBestg=Indexf(gBestpop,2);  
gBestfai=fai(gBestpop,:);  
gBestC=chrom(gBestpop,:,:);
gBest_x=chromx(gBestpop,:,:);
gBest_fit=fitness(gBestpop,:);
%----------------------------��ѭ����ʼ-----------------------------
for gen = 1:1:maxgen   
    for i = 1:1:vartotal 
        tmp=abs(chromx(1,gBestg,i)-oldchromx(1,gBestg,i));
        if tmp<1.0e-2
            tmp=1.0e-2;
        end
        max(i)=abs(fitness(1,gBestg)-oldfitness(1,gBestg))/tmp;
        for j = 1:1:popsize
            tmp=abs(chromx(j,gBestg,i)-oldchromx(j,gBestg,i)); 
            if tmp<1.0e-2
               tmp=1.0e-2;
            end
            if max(i)<abs(fitness(j,gBestg)-oldfitness(j,gBestg))/tmp
                max(i)=abs(fitness(j,gBestg)-oldfitness(j,gBestg))/tmp;
            end
            if min(i)>abs(fitness(j,gBestg)-oldfitness(j,gBestg))/tmp
                min(i)=abs(fitness(j,gBestg)-oldfitness(j,gBestg))/tmp;
            end
        end
    end
    %---------------ִ������λ��λ��ת���õ��µ���λ--------------------
    for i=1:1:popsize
        for j = 1:1:vartotal
            tmp=abs(chromx(i,gBestg,j)-oldchromx(i,gBestg,j));
            if tmp<1.0e-2
               tmp=1.0e-2;
            end
            grad=abs(fitness(i,gBestg)-oldfitness(i,gBestg))/tmp;
            tmp=abs(grad-min(j));
            if tmp<1.0e-2
               tmp=1.0e-2;
            end          
            rate(i,j)=tmp/abs(max(j)-min(j));            
            fai(i,j)=fai(i,j)+sign(chrom(i,1,j)*gBestC(1,2,j)-gBestC(1,1,j)*chrom(i,2,j))*(1-rate(i,j))*shiftstep*exp(-gen/maxgen);        
        end
    end
    %-----------------ִ������λ��λ����--------------------------
    Pm_rand = rand(popsize,vartotal);%������������������ʱȽϣ������Ƿ����
    for i=1:1:popsize
        for j=1:1:vartotal
            if (Pm(i)>Pm_rand(i,j))&&(i==gBestpop)
                fai(i,j)=0.5*pi-fai(i,j);
            end
        end
    end
    %-----------------���临��------������������������������Ӵ�Ⱦɫ��-------
    oldchrom=chrom;
    oldchromx=chromx;
    oldfitness=fitness;
    %---------------�����µ�����Ⱦɫ��-----------------
    chrom(:,1,:)=cos(fai(:,:));
    chrom(:,2,:)=sin(fai(:,:));
    %---------------��ռ�任-----------------------
    for i=1:1:2
        for j=1:1:vartotal
            chromx(:,i,j)=0.5*(var_range(j,2)*(1+chrom(:,i,j))+var_range(j,1)*(1-chrom(:,i,j)));
        end
    end
    %-----------������Ӧ��---��Ӧ�Ⱥ�����Shaffer's F6����-------------
    for i=1:1:popsize
        for j=1:1:2
            x1=chromx(i,j,1);
            x2=chromx(i,j,2);
            fitness(i,j)=0.5-((sin(sqrt(x1^2+x2^2)))^2-0.5)/(1+0.001*(x1^2+x2^2))^2;
        end
    end
    %----------------------������Ž⼰��Ӧ�Ա���----------------------
    [Bestf,Indexf]=sort(fitness,2);
    [BestF,IndexF]=sort(Bestf,1);
    Bestfit=BestF(popsize,2);  
    Bestpop=IndexF(popsize,2);
    Bestg=Indexf(Bestpop,2);  
    Bestfai=fai(Bestpop,:);  
    BestC=chrom(Bestpop,:,:); 
    Best_x=chromx(Bestpop,:,:);
    Best_fit=fitness(Bestpop,:);
    Badpop=IndexF(1,1); 
    %-----------------�����Ž��˻���ȡ���ϴ����Ž�------------------
    if Bestfit<gBestfit
        Bestfit=gBestfit;
        fai(Badpop,:)=gBestfai(1,:);
        chrom(Badpop,:,:)=gBestC(1,:,:);
        chromx(Badpop,:,:)=gBest_x(1,:,:);
        fitness(Badpop,:)=gBest_fit(1,:);
        gBestpop=Badpop;%���Ⱦɫ��ű�������
    end
    %---------------�����Ž���������Ž��滻------------------
    if Bestfit>=gBestfit
        gBestfit=Bestfit;  
        gBestpop=Bestpop; 
        gBestg=Bestg; 
        gBestfai=Bestfai;  
        gBestC=BestC; 
        gBest_x=Best_x;
        gBest_fit=Best_fit;
    end
    %---------------------��¼�Ż����---------------------
    result(gen)=gBestfit;
    iteration(gen)=gen;
    if result(gen)>0.995
        break;
    end
end
%-----------------��ѭ������-------------------
bestresult=result(gen);
iterationstep=iteration(gen);
bestresult
iterationstep
gBestg
figure(1)
plot(iteration,result)
