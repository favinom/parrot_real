close all
clear all

cd('/Users/favinom/projectsN/parrot_real/results');

folder = 'fracture-flow-3d-master/single/results/';

sub{1} = 'NCU_TW';

sub{2} = 'UIB/MPFA';
sub{3} = 'UIB/MVEM';
sub{4} = 'UIB/RT0';
sub{5} = 'UIB/TPFA';

sub{6} = 'UNICE-UNIGE/HFVCont';
sub{7} = 'UNICE-UNIGE/HFVDisc';
sub{8} = 'UNICE-UNIGE/VAGCont';
sub{9} = 'UNICE-UNIGE/VAGDisc';


sub{10} = 'USTUTT/MPFA';
%sub{11} = 'USI-UNIL/AD';
%sub{12} = 'USI-UNIL';
for ref=1:3
    figure(ref)
    
    for sub_i=1:length(sub)
        
        file = ['dot_refinement_',num2str(ref-1),'*'];
        
        orig=pwd;
        dove=[folder,sub{sub_i}];
        cd(dove);
        filename=dir(file);
        
        data{sub_i}=load(filename.name);
        
        cd(orig);
        
        
        
    end
    
    for ff=1:3
        subplot(3,2,2*ff-1);
        hold on
        for sub_i=1:length(sub)
           if (sub_i==12)
                plot(data{sub_i}(:,1),data{sub_i}(:,ff+1),'k-', 'LineWidth', 2)
           elseif (sub_i==11)
                plot(data{sub_i}(:,1),data{sub_i}(:,ff+1),'r-', 'LineWidth', 2)
           else
            plot(data{sub_i}(:,1),data{sub_i}(:,ff+1))
           end
        end
    end
    
end


for ref=1:3
    figure(ref)
    
    for sub_i=1:length(sub)
        
        file = ['dol_refinement_',num2str(ref-1),'*'];
        
        orig=pwd;
        dove=[folder,sub{sub_i}];
        cd(dove);
        filename=dir(file);
        
        data{sub_i}=load(filename.name);
        
        cd(orig);
        
    end
    
     for ff=1:3
         subplot(3,2,2*ff);
         hold on
         for sub_i=1:length(sub)
             if(sub_i==12)
              plot(data{sub_i}(:,2*ff-1),data{sub_i}(:,2*ff),'k-', 'LineWidth', 2)
             elseif (sub_i==11)
              plot(data{sub_i}(:,2*ff-1),data{sub_i}(:,2*ff),'r-', 'LineWidth', 2)
             else        
              plot(data{sub_i}(:,2*ff-1),data{sub_i}(:,2*ff))
             end
                          
         end
     end
     
end
