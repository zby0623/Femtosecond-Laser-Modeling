close all;
fs=19;
Nx=1000;Ny=10000;
T0=273.15;
dx=2E-5;
dy=2E-5;
T=(T0+37)*ones(Ny,Nx);
FWHM=40E-15;
sigma=FWHM/(2*sqrt(2*log(2)));
E=0.05E-3;
w0=2E-3;
P=E/(sigma*sqrt(2*pi));
dt=2E-15;
tfinal=300E-15;
h=20;
Tnew=T;
k=0.45;
rho=1050;
c=3421;
alpha=k/(rho*c);
t0=100E-15;
N_pulses=1;
Dt=tfinal/N_pulses;
t=0:dt:tfinal;
P_t=zeros(1,length(t));
for i=1:N_pulses
    P_t=P_t+P*exp(-(t-(t0+(i-1)*Dt)).^2/(2*sigma^2));
end
I=2*P_t/(pi*w0^2);
figure(1)
plot(t*1E15,I,'LineWidth',2.5);
xlabel('Time (fs)');
ylabel('Intensity (W/m^{2})');
set(gca,'Fontsize',fs);
set(gca,'FontWeight','Bold');
a=15;
f=10E-2;
T_f=zeros(round(tfinal/dt)+1,1);
T_f2=T_f;
T_1=T_f;
for t=0:dt:tfinal
    disp(t)
    for i=1:Ny-1
        for j=2:Nx-1
            w0_z=(f-i*dy)/f*w0;
            if i*dy==f
                w0_z=dy/f*w0;
            end
            q=I(round(t/dt)+1)*exp(-(j*dx-Nx/2*dx)^2/(w0_z^2))*exp(-a*i*dy);
            if i*dy>f
                q=0;
            end
            if i==1
                Tnew(i,j)=T(i,j)+dt*alpha*((T(i,j-1)+T(i,j+1)-2*T(i,j))/(dx^2)+(T(2,j)-T(1,j))/(dy^2))...
                    -h*(T(i,j)-T0)*dt/(rho*c)/dy+q*dt/(rho*c)/dy;
            else
                Tnew(i,j)=T(i,j)+dt*alpha*((T(i,j-1)+T(i,j+1)-2*T(i,j))/(dx^2)+(T(i-1,j)+T(i+1,j)-2*T(i,j))/(dy^2))...
                    +q*dt/(rho*c)/dy;
            end
        end
    end
    T=Tnew;
    T_f(round(t/dt)+1)=Tnew(round(f/dy),round(Nx/2));
    T_f2(round(t/dt)+1)=Tnew(round(f/2/dy),round(Nx/2));
    T_1(round(t/dt)+1)=Tnew(1,round(Nx/2));
end
%% plot
figure(2)
set(gcf,'Position',[100 100 700 500]);
imagesc((-Nx/2*dx+dx/2:dx:(Nx/2*dx-dx/2))*1E3,(dy/2:dy:Ny*dy-dy/2)*1E3,Tnew);
axis square;
axis tight;
xlabel('X (mm)');
ylabel('Z (mm)');
set(gca,'FontSize',fs);
set(gca,'FontWeight','Bold');
c=colorbar();
c.Label.String='Temperature (K)';
c.FontSize=21;
c.FontWeight='bold';
figure(3)
set(gcf,'Position',[100 100 600 500]);
plot((0:dt:tfinal)*1E15,T_1,'-o','LineWidth',2.5);
hold on;
plot((0:dt:tfinal)*1E15,T_f2,'-s','LineWidth',2.5);
hold on;
plot((0:dt:tfinal)*1E15,T_f,'-*','LineWidth',2.5);
hold off;
xlabel('Time (fs)');
ylabel('Temperature (K)');
legend('Z=0','Z=f/2', 'Z=f','Location','Northwest');
set(gca,'FontSize',fs);
set(gca,'FontWeight','Bold');
% axes('position',[.51 .3 .15 .35])
% box on % put box around new pair of axes
% t=0:dt:tfinal;
% Range_of_Interest=t>200E-15;
% plot(t(Range_of_Interest),T_1(Range_of_Interest),'LineWidth',2.5);
% hold on;
% plot(t(Range_of_Interest),T_f2(Range_of_Interest),'LineWidth',2.5);
% hold on;
% plot(t(Range_of_Interest),T_f(Range_of_Interest),'LineWidth',2.5);
% hold off;
% axis tight
% set(gca,'FontSize',16);
% set(gca,'FontWeight','bold');
% % axis([1.2 2.42 0.495 0.5])


