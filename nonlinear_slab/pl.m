figure(4001)
plot([0  1.7049    3.4099    5.1148    7.2784],[0 0.1 0.2 0.3 0.4],'b','LineWidth',2)
% xlim([minf,maxf])
xlabel('Lateral deflection','FontSize', 20)
ylabel('Load applied(N/mm^2)','FontSize', 20)
title('Load vs Deflection curve','FontSize', 20)
% legend('Col 600 X 600')
set(gcf,'PaperPosition',[0 0 10 7])
% print -djpeg -f4001 -r300

hold on;
plot([0 1.5172    3.0343    4.5515    6.0687    7.5858],[0 0.1 0.2 0.3 0.4 0.5],'r','LineWidth',2)
% xlim([minf,maxf])
xlabel('Lateral deflection','FontSize', 20)
ylabel('Load applied(N/mm^2)','FontSize', 20)
title('Load vs Deflection curve','FontSize', 20)
legend('Col 500 X 500','Col 600 X 600')
set(gcf,'PaperPosition',[0 0 10 7])
% print -djpeg -f4001 -r300

