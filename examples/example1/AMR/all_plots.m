close all

clear all

dol_0=csvread('1k/dol_refinement_0.csv');
dol_1=csvread('10k/dol_refinement_1.csv');
dol_2=csvread('100k/dol_refinement_2.csv');
k=0;
for (i=1:3)
    figure(i)
    plot(dol_0(:,i+k),dol_0(:,2*i) )
    k=k+1;
end
k=0;
for (i=1:3)
    figure(i)
    hold on
    plot(dol_1(:,i+k),dol_1(:,2*i) )
    k=k+1;
end
k=0;
for (i=1:3)
    figure(i)
    hold on
    plot(dol_2(:,i+k),dol_2(:,2*i) )
    k=k+1;
end

dot_0=csvread('1k/dot_refinement_0.csv');
dot_1=csvread('10k/dot_refinement_1.csv');
dot_2=csvread('100k/dot_refinement_2.csv');
k=0;
for (i=1:3)
    figure(i+3)
    plot(dot_0(:,1),dot_0(:,i+1) )
    k=k+1;
end
k=0;
for (i=1:3)
    figure(i+3)
    hold on
    plot(dot_1(:,1),dot_1(:,i+1) )
    k=k+1;
end
k=0;
for (i=1:3)
    figure(i+3)
    hold on
    plot(dot_2(:,1),dot_2(:,i+1) )
    k=k+1;
end