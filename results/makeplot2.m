close all
clear all

%cd('/Users/favinom/projectsN/parrot_real/results');

folder = 'fracture-flow-3d-master/regular/results/';

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
sub{11} = 'USI-UNIL';

for sub_i=1:length(sub)
    
    if ( isempty(sub{sub_i}) )
        continue;
    end
    
    
    orig=pwd;
    dove=[folder,sub{sub_i}];
    cd(dove);
    
    dirOut=dir('dot*');
    [c0,c1]=dirOut.name;
    file = c0;
    
    data{sub_i}=load(file);
    
    cd(orig);
    
end

for sub_i=1:length(data)
    if ( isempty(data{sub_i}) )
        continue;
    end
    
    sub_i
    size(data{sub_i})
    
    
end

for i=2:size(data{sub_i},2)-1
    
    figure(i-1);
    hold on
    for sub_i=1:length(data)
        if ( isempty(data{sub_i}) )
            continue;
        end
        
        if (sub_i==10)
            %aaaa=10;
            plot(data{sub_i}(:,1),data{sub_i}(:,i),'*')
        elseif (sub_i==11)
             plot(data{sub_i}(:,1),data{sub_i}(:,i),'^')
        elseif (sub_i==12)
            plot(data{sub_i}(:,1),data{sub_i}(:,i),'o')
        else
            plot(data{sub_i}(:,1),data{sub_i}(:,i))
        end
    end
end
