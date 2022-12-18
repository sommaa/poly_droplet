%                              ________  ________  ___           ___    ___ ________                                   %
%                             |\   __  \|\   __  \|\  \         |\  \  /  /|\   ___ \                                  %
%                             \ \  \|\  \ \  \|\  \ \  \        \ \  \/  / | \  \_|\ \                                 %
%                              \ \   ____\ \  \\\  \ \  \        \ \    / / \ \  \ \\ \                                %
%                               \ \  \___|\ \  \\\  \ \  \____    \/  /  /   \ \  \_\\ \                               %
%                                \ \__\    \ \_______\ \_______\__/  / /      \ \_______\                              %
%                                 \|__|     \|_______|\|_______|\___/ /        \|_______|                              %
%                                                              \|___|/                                                 %
%                                            Author: Andrea Somma;                                                     % 
%                                            Politecnico of Milan 2021-2022                                            % 
%                                                                                                                      %

clc; clear all

%-----------------------------------------------------------------------------------------------------------------------
%% initial data
%-----------------------------------------------------------------------------------------------------------------------
D = 4e-9;
Cw0 = 55*1e-3; %mol/m3;
r0 = 1e-3; %m
K = 0.595;
Cm = 25.*1e-3;
t00 = 1/K;
t1t = Cm/Cw0;
delta0 = 0;
Y = (1/t1t + 1 - K)^(1/3) -1;
rdF = Y*r0 + r0;

L = 1; %adim
tf = 3; %adim

%-----------------------------------------------------------------------------------------------------------------------
%% preprocessing
%-----------------------------------------------------------------------------------------------------------------------
nsteps = 15;
h = L/nsteps;
M = speye(nsteps+1,nsteps+1);
M(1,1) = 0;

options = odeset('Mass',M,'MaxStep',1e-5,'RelTol',1e-4,'AbsTol',1e-6,'InitialStep',1e-15); %NON CAMBIARE


init = [linspace(t00,t1t,nsteps),1e-8]'; % se lo cambi ti penti
[t,X] = ode15s(@(t,X) ODESYS(t,X,nsteps,h,t1t,Y),0:0.01:tf,init,options);

%-----------------------------------------------------------------------------------------------------------------------
%% plot 1D
%-----------------------------------------------------------------------------------------------------------------------
close all
A = 1:1:150; B = 201:2:499; C = 499:4:size(X,1);
plotter = [A,B,C];

for k=2:2:size(X,1)
    hold off
    plot([0,0.001005],[X(k,1).*Cw0 X(k,1).*Cw0], 'linewidth',2)
    hold on;
    plot((h:h:L).*X(k,end).*(rdF-r0)+r0, X(k,1:end-1).*Cw0, 'linewidth',2); 
    ylim([0.02 0.06])
    xlim([0 X(end,end).*(rdF-r0)+r0*1.1])
    xlabel("radius [m]")
    ylabel("Cw [mol/m3]")  
    message = sprintf('time=%g s', round(t(k),4).*(rdF-r0).^2./D);
    time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','k');
    drawnow
    delete(time)
end

%-----------------------------------------------------------------------------------------------------------------------
%% plot radius and concentration
%-----------------------------------------------------------------------------------------------------------------------
figure(2)
yyaxis left
set(gcf,'position',[50,50,900,400])
plot(t.*(rdF-r0).^2./D,X(:,end).*(rdF-r0)+r0)
ylabel("radius [m]")
hold on

yyaxis right
plot(t.*(rdF-r0).^2./D,X(:,1).*Cw0)
hold on
title("radius and internal water concentration variation in time")
ylabel("Cw [mol/m3]")
xlabel("time [s]")
message = sprintf('Diff=%g m2/s, r0=%g m,\nCw0=%g mol/m3, K=%g,\n Cm=%g mol/m3', [D r0 Cw0 K Cm]);
time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','k','BackgroundColor','w');

%-----------------------------------------------------------------------------------------------------------------------
%% plot 2D
%-----------------------------------------------------------------------------------------------------------------------
video_name = 'polydrop.mp4';
videompg4 = VideoWriter(video_name, 'Motion JPEG AVI');
open(videompg4);

figure(3)
set(gcf,'position',[1000,60,700,600])
set(gca,'XColor', 'none','YColor','none');

% interp grid
x = linspace(-r0*1.5, r0*1.5, 250);
y = linspace(-r0*1.5, r0*1.5, 250);
R = sqrt(x.^2 + (y.').^2);

for i=1:numel(plotter)
    % 2d interpolation
    k = plotter(i);
    hold off
    rad1D = [0 0.001,(h:h:L).*X(k,end).*(rdF-r0)+r0];
    C = [X(k,1).*Cw0 X(k,1).*Cw0, X(k,1:end-1).*Cw0];
    z = interp1(rad1D',C',R);
    % time
    message = sprintf('time=%g s \nrd=%g mm\nCw=%g mol/m3', round(t(k).*(rdF-r0).^2./D,2),round((X(k,end) ...
        .*(rdF-r0)+r0)*1e3,3),round(X(k,1).*Cw0,4));
    time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','k','BackgroundColor' ...
        ,'w','Interpreter','tex','LineWidth',0.1);
    %plot C
    imagesc(z,"Interpolation","bilinear");
    % visual stuff
    colorbar
    clim([Cm Cw0])
    colormap(jet(8e2))
    title("Water concentration mol/m3")
    axis off
    drawnow
    hold on

    frame = getframe(gcf);
    writeVideo(videompg4,frame);
    delete(time)

end

close(videompg4);

%-----------------------------------------------------------------------------------------------------------------------
%% function
%-----------------------------------------------------------------------------------------------------------------------
function dxdt = ODESYS(~,X,nsteps,h,t1t,Y)
    
    dxdt = zeros(1,nsteps+1);
    
    
    % concentration vector
    theta = X(1:nsteps);
    % int deltaquadro
    dxdt(nsteps+1) = - (theta(nsteps)-theta(nsteps-1))/h*2/t1t;
    %deltaquadro(t)
    deltaquadro = X(nsteps+1);
    
    delta = sqrt(deltaquadro);

    %integrazione int thetaw trapezi
    int = 0;
    for i = 1:nsteps-1
        int = int + ((h*i*delta*Y+1)^2*theta(i) + (h*(i+1)*delta*Y+1)^2*theta(i+1))*h/2;
    end
    
    %thetaw
    thetaw = 1 - 3*delta*Y*int;
    
    % condizione west
    dxdt(1) = thetaw - theta(1);
    
    % d(1/delta/Y)/dR neglected
    for i = 2:nsteps-1
        dxdt(i) = ( 2/(i*h+1/delta/Y)*(theta(i+1)-theta(i-1))/h + (theta(i+1)-2*theta(i)+theta(i-1))/h^2 + ...
            i/2*dxdt(nsteps+1)*(theta(i+1)-theta(i-1))/2 )/deltaquadro;
    end

    % condizione east
    dxdt(nsteps) = 0;
    dxdt = dxdt';
end
