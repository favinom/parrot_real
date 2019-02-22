function A=moveNodesPerfect(A)

thickness=0.5e-4;
q=(0.5+0.75)/2;
centers=[0.5 q 0.75];

P=A(:,1);
points=unique(sort(P));
numEl=(length(points)-4)/2;

x1=linspace(0,centers(1)-thickness,numEl+1);
x2=linspace(centers(1)+thickness,centers(2)-thickness,numEl/4+1);
x3=linspace(centers(2)+thickness,centers(3)-thickness,numEl/4+1);
x4=linspace(centers(3)+thickness,1,numEl/2+1);

x=[x1 x2 x3 x4];
    
    for coord=1:3
        for i=1:length(P)
            
            A(i,coord)=x(points==A(i,coord));
        
         end
    end
