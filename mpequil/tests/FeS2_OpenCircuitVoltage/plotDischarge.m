% home = cd;
% cd /home/clueth/Matlab/datagrabber/; h = datagrabber('LiSi_Discharge_Sharma1977',-1)
% cd(home)
h = [];
%title{1} = 'Experiments of Sharma et al. 1977';
title{2} = 'Cantera Simulation';
files = dir('*ocv*');
figure(1)
% hold on
n = length(files);
colors = hsv(n);
%current = [0.2455 0.491 0.982 1.964]
for i = 1:n
    [nFe, nS, nLi, capacity, voltage, junk] = textread(files(i).name,'','delimiter',',','headerlines',3);
    xLi = nLi./(nFe);
    %     voltage = 0.194+voltage-voltage(1);
    %voltage = 0.04*current(i)+voltage;
    hh(i) = plot(xLi,voltage,'k-*');
    %    hh(i) = plot(capacity/1.825,voltage,'--','Color',colors(i,:));
    hold on
end
%set(2,'PaperSize',[8.5 11])
%xlabel('Capacity [Ah/g Li_{13}Si_4]','FontSize',20)
%xlabel('Capacity [Amp-h]','FontSize',20)
xlabel('l: Li_lFeS2','FontSize',20)
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
% if isempty(h)
%     legend(hh(1),title{2},'FontSize',20,'location','NorthWest')
% else
%     legend([h(1) hh(1)],title{1},title{2},'FontSize',20,'location','NorthWest')
% end

%axis([0 3.5 0 0.6])
hold off

print(1,'-djpeg99','FeS2_Discharge.jpg')