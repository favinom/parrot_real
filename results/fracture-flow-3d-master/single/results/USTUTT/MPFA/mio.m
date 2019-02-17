
clear all

syms x;

phi{1}=1-x;
phi{2}=x;


for i=1:2
    for j=1:2
        
        temp=diff(phi{i});
        
        l(i,j)=-int(phi{j}*temp,0,1);
        
    end
end