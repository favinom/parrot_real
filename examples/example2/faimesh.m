clear all
close all

fid=fopen('mesh0_reg_0001_mesh.xda');
wfid=fopen('adapt1.xda','w');

dummy=textscan(fid,'%s',1,'delimiter','\n');
fprintf(wfid,[dummy{1}{1},'\n']);
elemStr=textscan(fid,'%s',1,'delimiter','\n');
fprintf(wfid,[elemStr{1}{1},'\n']);

C = strsplit(elemStr{1}{1},' ');
elemNum=str2num(C{1});

nodeStr=textscan(fid,'%s',1,'delimiter','\n');
fprintf(wfid,[nodeStr{1}{1},'\n']);

C = strsplit(nodeStr{1}{1},' ');
nudeNum=str2num(C{1});


for i=1:14
    InputText=textscan(fid,'%s',1,'delimiter','\n');
    fprintf(wfid,[InputText{1}{1},'\n']);
end

for i=1:elemNum
    InputText=textscan(fid,'%s',1,'delimiter','\n');
    fprintf(wfid,[InputText{1}{1},'\n']);
end

P=zeros(nudeNum,3);

for i=1:nudeNum
    InputText=textscan(fid,'%s',1,'delimiter','\n');
    C = strsplit(InputText{1}{1},' ');
    P(i,1)=str2double(C{1});
    P(i,2)=str2double(C{2});
    P(i,3)=str2double(C{3});
end

A=moveNodes(P);

for i=1:nudeNum
    fprintf(wfid,[num2str( A(i,1) ),' ']);
    fprintf(wfid,[num2str( A(i,2) ),' ']);
    fprintf(wfid,[num2str( A(i,3) ),'\n']);
end

while (~feof(fid))
    InputText=textscan(fid,'%s',1,'delimiter','\n');
    fprintf(wfid,[InputText{1}{1},'\n']);
end


fclose(wfid);
return

plot3(A(:,1),A(:,2),A(:,3),'*');

return



