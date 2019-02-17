%script to compute parameters
function[K_eq,k_eq]= compute_param(e,K,k)
K_eq=zeros(length(K),1);
k_eq=zeros(length(k),1);
a=zeros(length(k),1);
e_tot=zeros(length(e)+1,1);
for i=1:length(e)
    e_tot(i)=e(end+1-i);
end
e_tot(length(e)+1)=1.0;
dim_f=2;
for i=1:length(K)
    K_eq(i)= K(i)/e(i);
    a(i)=nthroot(e(i),3-dim_f);
    a_v=a(i)
    e_v=e_tot(dim_f+2)
    k_eq(i)=k(i)*a(i)/(2*e_tot(dim_f+2));
    k_v=k_eq(i)
    dim_f=dim_f-1;
end
end
