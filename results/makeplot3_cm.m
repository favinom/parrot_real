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
sub{11} = 'USI-UNIL';

for sub_i=1:length(sub)
    
    if ( isempty(sub{sub_i}) )
        continue;
    end
    
    
    orig=pwd;
    dove=[folder,sub{sub_i}];
    cd(dove);
    
    dirOut=dir('dot*0*');
    
    if ( length(dirOut) ~= 0)
        filename=dirOut.name;
        
        if (sub_i==11)
            data{sub_i}=csvread(filename,2);
            data{sub_i}=[data{sub_i}(:,1) data{sub_i}(:,2:9)./data{sub_i}(:,10:end)];
        else
            data{sub_i}=load(filename);
        end
    end
    
    cd(orig);
    
end


for sub_i=1:length(data)
    if ( isempty(data{sub_i}) )
        continue;
    end
    
    %sub_i
    %size(data{sub_i})
    
    
end


for fig=2:size(data{sub_i},2)
    figure(fig-1)
    hold on
    counter=0;
    for sub_i=1:length(data)
        if ( isempty(data{sub_i}) )
            continue;
        end
        
        if (sub_i==11)
            plot(data{sub_i}(:,1),data{sub_i}(:,fig),'*')
        else
            plot(data{sub_i}(:,1),data{sub_i}(:,fig))
        end
        counter=counter+1;
        l{counter}=sub{sub_i};
    end
    
    for i=1:counter
        legend(l)
    end
    
end