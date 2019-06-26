function A=moveNodes(A)

q=(0.5+0.75)/2;
centers=[0.5 0.5 0.5;
         0.75 0.75 0.75;
         q q q];
     
thickness=0.0001;

[mC,~]=size(centers);

for centIter=1:mC
    
    for coord=1:3
        
        P=A(:,coord);
        S=sort(P);
        
        whereP=find(S>centers(centIter,coord)+thickness/2);
        valueP=S(whereP(1));
        
        whereN=find(S<centers(centIter,coord)-thickness/2);
        valueN=S(whereN(end));
        
        
        for i=1:length(P)
            if (abs(valueP-P(i))<1e-8)
                P(i)=centers(centIter,coord)+thickness/2;
            end
            if (abs(valueN-P(i))<1e-8)
                P(i)=centers(centIter,coord)-thickness/2;
            end
        end
        
        A(:,coord)=P;
        
    end
end

%P=A(:,1);
%unique(sort(P))