%script to compute parameters
function[K_eq,k_eq]= compute_param(e,K,k,dim)
K_eq=zeros(length(K),1);
k_eq=zeros(length(k),1);
a=zeros(length(k),1);
e_tot=zeros(length(e)+1);
e_tot=e;
e_tot(length(e)+1)=1.0;
for i=1:length(K)
    K_eq(i)= K(i)/e_tot(i);
    if(dim-1>-1)
        a(i)=nthroot(e_tot(i),dim-1);
    end
    k_eq(i)=k(i)*a(i)/(2*e_tot(i+1));
end
end
