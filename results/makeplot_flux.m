close all
clear all

%cd('/Users/favinom/projectsN/parrot_real/results');

folder = 'fracture-flow-3d-master/small_features/results/';

sub{1} = [];

sub{2} = 'UIB/MPFA';
sub{3} = 'UIB/MVEM';
sub{4} = 'UIB/RT0';
sub{5} = 'UIB/TPFA';

sub{6} = 'UNICE-UNIGE/HFVCont';
sub{7} = 'UNICE-UNIGE/HFVDisc';
sub{8} = 'UNICE-UNIGE/VAGCont';
sub{9} = 'UNICE-UNIGE/VAGDisc';


sub{10} = 'USTUTT/MPFA';
%sub{11} = 'USI-UNIL';

for sub_i=1:length(sub)
    
    if ( isempty(sub{sub_i}) )
        continue;
    end
    
    
    orig=pwd;
    dove=[folder,sub{sub_i}];
    cd(dove);
    
    dirOut=dir('resul*');
    
    if ( length(dirOut) ~= 0)
        filename=dirOut.name;
        
        data{sub_i}=load(filename);
        
    end
    
    cd(orig);
    
end

l=[];
for sub_i=1:length(data)
    if ( isempty(data{sub_i}) )
        continue;
    end
    
    l=[l data{sub_i}(2,end)];
    %size(data{sub_i})
    
    
end
