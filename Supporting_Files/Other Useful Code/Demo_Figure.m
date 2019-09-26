set(0,'defaultAxesFontSize',20);
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 2)
close all
clear all

% seed = rng, save seed
% load seed
% rng(seed)

%set up signals, filters etc
dt = 0.0006;
t = 0:dt:1023*dt;
N = length(t);
fs = 1/dt;
fmin = 1/(t(end)-t(1));
T = 0.105;
w = 2*pi/T;

%generate signals
phi0 = -2.3;
y0 = 0.4*(cos(w*t-phi0) + 0.003*rand(1,N));
y0 = y0 + 0.5;
y1 = 0.4*(cos(w*t-phi0-pi/3) + 0.003*rand(1,N));
y1 = y1 + 0.5;

%digitise signals
nDigi = 10;
y0 = round(y0*(2^nDigi-1))/(2^nDigi-1);
y1 = floor(y1*(2^nDigi-1))/(2^nDigi-1);

%misc parameters
W = hanning(N,'periodic')';
% W = W.^0;
CN = 2;%correction_factor(W,length(W));
k = 0:N-1;
F = linspace(0,fs*(N-1)/N,N);
dF = diff(F); dF = dF(1);
fLow = 0;
fHi = 12;

%time domain windowing
Y0a = fft(y0.*W);
Y1a = fft(y1.*W);
phi0a = wrapToPi(angle(Y0a));
phi1a = wrapToPi(angle(Y1a));
R0a = 10*abs(Y0a)/N;
R1a = 10*abs(Y1a)/N;

%find f0
[~,I] = max(abs(Y0a(3:round(N/2)))); %ignoring DC
I = I + 2; %correcting for omission of DC
kEst0a = I - 1 + CN*real((Y0a(I-1)-Y0a(I+1))/(2*Y0a(I)-Y0a(I-1)-Y0a(I+1)));
[~,I] = max(abs(Y1a(3:round(N/2))));
I = I + 2;
kEst1a = I - 1 + CN*real((Y1a(I-1)-Y1a(I+1))/(2*Y1a(I)-Y1a(I-1)-Y1a(I+1)));
kEsta = (kEst0a+kEst1a)/2;

%find dphi0a
dphi0a = interp1(1:N,wrapToPi(phi1a-phi0a),kEsta,'linear');

%freq domain windowing
Y0b = sdft_plug(y0);
Y1b = sdft_plug(y1);
phi0b = angle(Y0b);
R0b = 10*abs(Y0b)/N;
phi1b = angle(Y1b);
R1b = 10*abs(Y1b)/N;
dummyA = [NaN Y0b(1:end-1)];
dummyB = [Y0b(2:end) NaN];
Y0c = -0.25*dummyA + 0.5*Y0b - 0.25*dummyB;
R0c = 10*abs(Y0c)/N;
phi0c = angle(Y0c);
dummyA = [NaN Y1b(1:end-1)];
dummyB = [Y1b(2:end) NaN];
Y1c = -0.25*dummyA + 0.5*Y1b - 0.25*dummyB;
R1c = 10*abs(Y1c)/N;
phi1c = angle(Y1c);

%find f0
[~,I] = max(abs(Y0c(3:round(N/2))));
I = I + 2;
kEst0c = I - 1 + CN*real((Y0c(I-1)-Y0c(I+1))/(2*Y0c(I)-Y0c(I-1)-Y0c(I+1)));
[~,I] = max(abs(Y1c(3:round(N/2))));
I = I + 2;
kEst1c = I - 1 + CN*real((Y1c(I-1)-Y1c(I+1))/(2*Y1c(I)-Y1c(I-1)-Y1c(I+1)));
kEstc = (kEst0c+kEst1c)/2;

if abs(kEst1c-kEst0c)/kEst0c > 0.01
    disp('Help');
end

%find dphi0c
dphi0c = interp1(1:N,wrapToPi(phi1c-phi0c),kEstc,'linear');


%time domain windowing

figure('units','normalized','outerposition',[0 0 0.5 0.5]);
pause(0.5);

axes('position',[0.1 0.85 0.2 0.11]);
plot(t,y0,t,y1);
% legend('0','1');
xlabel('Time (s)');
set(gca,'xtick',[0 0.2 0.4 0.6]);
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('Signal (a.u.)');
yyaxis right
set(gca, 'YColor', 'k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.005,0.875,'(a)','fontsize',20);

axes('position',[0.1 0.570 0.2 0.11]);
plot(t,W.*y0,t,W.*y1);
set(gca,'xtick',[0 0.2 0.4 0.6]);
% legend('0','1');
xlabel('Time (s)');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('Signal (a.u.)');
yyaxis right
set(gca, 'YColor', 'k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.005,0.875,'(b)','fontsize',20);

axes('position',[0.1 0.30 0.2 0.11]);
plot(k(3:end-1),R0a(3:end-1),'-',k(3:end-1),R1a(3:end-1),'o','markersize',8);
line([kEsta kEsta],[0 max(R0a(3:end))],'linestyle','-','color','k','linewidth',1);
set(gca,'xlim',[fLow fHi]);
set(gca,'xticklabel',[]);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\fontname{Times}R} (a.u.)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,1.1,'(c)','fontsize',20);

axes('position',[0.1 0.175 0.2 0.11]);
plot(k(3:end-1),phi0a(3:end-1),'-o',k(3:end-1),phi1a(3:end-1),'-o','markersize',8);
line([kEsta kEsta],get(gca,'ylim'),'linestyle','-','color','k','linewidth',1);
set(gca,'xlim',[fLow fHi]);
set(gca,'xticklabel',[]);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\phi} (rad)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,2.3,'(d)','fontsize',20);

axes('position',[0.1 0.05 0.2 0.11]);
plot(k(3:end-1),wrapToPi(phi1a(3:end-1)-phi0a(3:end-1)),'-o','markersize',8);
line([kEsta kEsta],get(gca,'ylim'),'linestyle','-','color','k','linewidth',1);
xlabel('$$k$$','Interpreter','Latex');
set(gca,'xlim',[fLow fHi]);
set(gca,'xtick',0:3:12);
set(gca,'ylim',[-1.2 -0.9]);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\Delta}{\it{\phi}} (rad)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,-0.935,'(e)','fontsize',20);

%freq domain filtering

axes('position',[0.4 0.85 0.2 0.11]);
plot(t,y0,t,y1);
set(gca,'xtick',[0 0.2 0.4 0.6]);
% legend('0','1');
xlabel('Time (s)');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('Signal (a.u.)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.005,0.875,'(f)','fontsize',20);


axes('position',[0.4 0.635 0.2 0.11]);
plot(k(3:end-1),R0b(3:end-1),'-',k(3:end-1),R1b(3:end-1),'o','markersize',8);
set(gca,'xticklabel',[]);
set(gca,'xlim',[fLow fHi]);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\fontname{Times}R} (a.u.)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,1.75,'(g)','fontsize',20);

axes('position',[0.4 0.51 0.2 0.11]);
plot(k(3:end-1),phi0b(3:end-1),'-o',k(3:end-1),phi1b(3:end-1),'-o','markersize',8);
xlabel('$$k$$','Interpreter','Latex');
set(gca,'xlim',[fLow fHi]);
set(gca,'xtick',0:3:12);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\phi} (rad)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,1.2,'(h)','fontsize',20);

axes('position',[0.4 0.30 0.2 0.11]);
plot(k(3:end-1),R0c(3:end-1),'-',k(3:end-1),R1c(3:end-1),'o','markersize',8);
line([kEstc kEstc],[0 max(R0c)],'linestyle','-','Color','k','linewidth',1);
set(gca,'xlim',[fLow fHi]);
set(gca,'xticklabel',[]);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\fontname{Times}R} (a.u.)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,1.1,'(i)','fontsize',20);

axes('position',[0.4 0.175 0.2 0.11]);
plot(k(3:end-1),phi0c(3:end-1),'-o',k(3:end-1),phi1c(3:end-1),'-o','markersize',8);
line([kEstc kEstc],get(gca,'ylim'),'linestyle','-','color','k','linewidth',1);
set(gca,'xlim',[fLow fHi]);
set(gca,'xticklabel',[]);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\phi} (rad)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,2.3,'(j)','fontsize',20);

axes('position',[0.4 0.05 0.2 0.11]);
plot(k(3:end-1),wrapToPi(phi1c(3:end-1)-phi0c(3:end-1)),'-o','markersize',8);
line([kEstc kEstc],get(gca,'ylim'),'linestyle','-','Color','k','linewidth',1);
set(gca,'xlim',[fLow fHi]);
set(gca,'ylim',[-1.2 -0.9]);
xlabel('$$k$$','Interpreter','Latex');
set(gca,'xtick',0:3:12);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\Delta}{\it{\phi}} (rad)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,-0.935,'(k)','fontsize',20);

annotation(gcf,'arrow',[0.125 0.125],[0.837 0.694]);
annotation(gcf,'arrow',[0.125 0.125],[0.837 0.694]-0.275);
annotation(gcf,'arrow',[0.4275 0.4275],[0.835 0.76]);
annotation(gcf,'arrow',[0.435 0.435],[0.835 0.76]-0.3375);
annotation(gcf,'textbox',...
    [0.076 0.736 0.048 0.0582],...
    'String',{'apply','window'},...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'HorizontalAlignment','right');

annotation(gcf,'textbox',...
    [0.095 0.478 0.030 0.032],...
    'String',{'DFT'},...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.399 0.785 0.0283 0.032],...
    'String',{'DFT'},...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right');

annotation(gcf,'textbox',...
        [0.3975 0.455 0.035 0.032],...
    'String',{'apply','filter'},...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right');

dummyString = ['$$\hat{k}$$ = ' num2str(kEsta,'%4.3f')];
annotation(gcf,'textbox',...
        [0.185 0.373 0.113 0.032],...
    'String',dummyString,...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right',...
    'Interpreter','Latex');

dummyString = ['$$\hat{k}$$ = ' num2str(kEstc,'%4.3f')];
annotation(gcf,'textbox',...
        [0.485 0.373 0.113 0.032],...
    'String',dummyString,...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right',...
    'Interpreter','Latex');

dummyString = ['$$\Delta{\hat{\phi}}$$ = ' num2str(dphi0a,'%4.3f')];
annotation(gcf,'textbox',...
        [0.185 0.124 0.113 0.032],...
    'String',dummyString,...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right',...
    'Interpreter','Latex');

dummyString = ['$$\Delta{\hat{\phi}}$$ = ' num2str(dphi0c,'%4.3f')];
annotation(gcf,'textbox',...
        [0.485 0.124 0.113 0.032],...
    'String',dummyString,...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right',...
    'Interpreter','Latex');


set(gcf,'color','w')

disp(['f0a = ' num2str(dF*kEsta)]);
disp(['f0c = ' num2str(dF*kEstc)]);
disp(['% error on dphi = ' num2str(abs(pi/3+dphi0c)/(pi/3))]);
