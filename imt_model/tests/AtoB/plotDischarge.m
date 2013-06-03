home = cd;
h = [];
% cd /home/clueth/Matlab/datagrabber/; h = datagrabber('LiSi_Discharge_Sharma1977',-1)
cd /home/clueth/Matlab/datagrabber/; h = datagrabber('fes2_discharge',-1)
set(h,'Color','k')
cd(home)
title{1} = 'Open Circuit Voltage';
title{2} = 'Cantera, single reaction';
title{3} = 'Cantera, multiple reactions';
files = dir('*Hist*');
figure(2)
hold on
n = length(files);
colors = hsv(n);
current = [0.2455 0.491 0.982 1.964]
for i = 1:n
    [time, charge, capacity, voltage] = textread(files(i).name,'','delimiter',',','headerlines',5);
    %     voltage = 0.194+voltage-voltage(1);
    %voltage = 0.04*current(i)+voltage;
    hh(i) = plot(-capacity/1.8766/2*100,voltage-2+1.908,'--','Color',colors(i,:));
    %    hh(i) = plot(capacity/1.825,voltage,'--','Color',colors(i,:));
    hold on
end
%set(2,'PaperSize',[8.5 11])
%xlabel('Capacity [Ah/g Li_{13}Si_4]','FontSize',20)
%xlabel('Capacity [Amp-h]','FontSize',20)
xlabel('% Utilization FeS_2','FontSize',20)
ylabel('Volts vs. Li(l)','FontSize',20)
set(gca,'FontSize',16)
% %text(0.52,0.17,'i = 0.13 A/g Li_{13}Si_4','FontSize',16,'Color',colors(1,:))
% %text(0.15,0.38,'i = 1.08 A/g Li_{13}Si_4','FontSize',16,'Color',colors(n,:))
% text(1,0.17,'i = 0.13 Amps','FontSize',16,'Color',colors(1,:))
% text(0.5,0.38,'i = 1.08 Amps','FontSize',16,'Color',colors(n,:))
% %text(1.4,0.07,'T = 402 ^oC','FontSize',16,'Color','k')
% text(2.4,0.12,'1.825 g Li_{13}Si_4','FontSize',16,'Color','k')
% text(2.5,0.07,'T = 402 ^oC','FontSize',16,'Color','k')

% 
if isempty(h)
    legend(hh(1),title{2},'FontSize',20,'location','SouthWest')
else
    legend([h(1) hh(1) hh(2)],title{1},title{2},title{3},'FontSize',20,'location','SouthWest')
end

%axis([-0.35 0 1.975 1.995])
axis([0 65 0 1.92])
hold off

%g = get(gca,'xlabel');set(g,'FontSize',20)
%set(2,'paperpositionmode','auto'); 
%print(2,'-djpeg99','FeS2_Discharge.jpg')
%print(2,'-dpng','FeS2_Discharge.png')
print(2,'-depsc','FeS2_Discharge.eps')
