%one dimensional model of main stem that evolves through time according to 
% stream power model of bedrock incision  dz/dt = U - K*A^m*S^n over the domain Xcrit < x < L
% (the fluvial part of a longitudinal profile) and simple prescribed uplift fields.  
%No accounting for sediment flux and its role 
% erosion and sedimentation on the bed is included.

% Currently a very simple rule is used to turn on waterfall processes or other 
% sediment transport related processes that would systematically perturb the erosion 
% rate above or below that calculated using stream power.

%uses function from matlab file exchange:
%       cbrewer()
%Code by mostly Scott McCoy, figures by Sophie Rothman


close all
clear all

%% set parameter values
% model time
tf = 15e6;              %final run time in Myr
dt = 10;                %time step (yr)
ta = 0:dt:tf;

% define model grids 
L = 30e3;              % river length [m]
zb = 0;                 %base level of mainstem [m]
dx = 50;                %grid spacing
xc = 400;               %set critical distance from ridge to define rivers (ie where fluvial processes turn on [m]
x_d = (xc:dx:L);     %build x array of distance from divide
x = 0:dx:L-xc;        %%build x array of distance from outlet or baselevel

% set parameter values
%stream power erosion
ka = 6.69;          %hack coefficient
h = 1.8;            %inverse hack exponent
K = 5e-7;           %erodibility
n = 1;              %slope exponent
m = 0.45;           %drainage area exponent
theta = m/n;        % channel concavity

%uplift parameters
U = 80e-6;          %uplift rate [m/yr]
%Ufield = linspace(U,0,length(x)); %change the uplift rate to a tilt
Ufield = ones(1,length(x))*U;      %uniform U
%Ufield= ones(1,length(x))*U*3;      %step increase in U
Ufield(end) = 0;                   %set U =0 at fixed baselevel
counter=1;
%waterfalls
waterfalls_ms = 1  ;           %turn on autogenic waterfalls in mainstem 
Fms = 2  ;                  % factor to change efficientcy of stream power rate
Sc = 0.05    ;                 %Slope above which waterfalls are turned on

if waterfalls_ms
   xcrit=(((U/(K*(ka^m)))^(1/(h*m)))*((1/Sc)^(n/(h*m))));
   zcrit= (((U/K)^(1/n))*((ka^(-theta))/(1-(h*theta)))*(((L)^(1-(h*theta))) - (x_d((int16(((xcrit-xc)/dx))-1))^(1-(h*theta)))));   
end


%plots
plotevery = 50000;      %plot every so many yr
saveevery = .2e6;       %save frames every so many yr
zmax = 3000;             %max elevation to plot [m]

blues = [16,78,139; 30,144,255; 0,191,255]./255;  
%colors = cbrewer('seq', 'Blues', round(tf/saveevery)+10);   %color map for ploting multiple times on one plot


%initialize other  vars
count = 9;   %counter to step through color array start at 9 to skip the first four that are too light to see
t = 0; %initial time
iteration = 0;
wfix =[]
wfixTVD =[]
Nnodes = length(x);
dzdx = zeros(1,Nnodes);
dzdx0=zeros(1,Nnodes);
dzdx2=zeros(1,Nnodes);
dzdx_wf=zeros(1,Nnodes);
zwf=zeros(1,Nnodes);
counter_v=zeros(1,Nnodes);
counter_d=zeros(1,Nnodes);
savedzdx_v=zeros(int8(tf/saveevery),Nnodes+1);
savedzdx_d=zeros(int8(tf/saveevery),Nnodes);
savez_v=zeros(int8(tf/saveevery),Nnodes+1);
savez_d=zeros(int8(tf/saveevery),Nnodes);
savewf_v=zeros(int8(tf/saveevery),Nnodes+1);
savewf_d=zeros(int8(tf/saveevery),Nnodes);

s1 = zeros(1,Nnodes);  %
start_waterfalls=1;



%% set initial conditions MAKE SURE ONLY ONE IS UNCOMMENTED

% %% ---- analytical steady state initial condition
%z0 = ((U/K)^(1/n))*((ka^-theta)/(1-h*theta)).*((x_d(end)^(1-h*theta)) - (x_d.^(1-h*theta)));
    %z0=flip(z0);    %flip z0 such that  x(1) is the outlet instead of the divide

% ---- Plateau (constant elevation profile)
% z_plat = 0;	       
% z0 = ones(size(x)).*z_plat;
% z0(end) = 0;

% %% ----inclined plane
% islope = 0.05
% z0 = (L + (-1.*x_d)).*islope;


% %% ----instantaneous uniform uplift of steady state profile
Upulse = 200; %200;        %amount to uplift in instantaneous pulse [m] set to 0 if no pulse
Pulse_interval=.5e6;
pulse_per_unit=(Upulse/Pulse_interval)*dt;
Pulseovertime=zeros(length(ta), 1);
Pulseovertime(1:(Pulse_interval/dt), 1)=ones((Pulse_interval/dt), 1)*pulse_per_unit;



z0 = ((U/K)^(1/n))*((ka^-theta)/(1-h*theta)).*((x_d(end)^(1-h*theta)) - (x_d.^(1-h*theta)));
% z=z0

% % %% ----instantaneous tilt uplift of steady state profile
% Upulse =  500;        %amount to uplift in instantaneous pulse [m] set to 0 if no pulse
% z0 = ((U/K)^(1/n))*((ka^-theta)/(1-h*theta)).*((x_d(end)^(1-h*theta)) - (x_d.^(1-h*theta)));
% z0(1:end-1) = z0(1:end-1) + flip(linspace(0,Upulse,Nnodes-1));

%Waterfall Steady State
start_waterfalls=1;
if waterfalls_ms==1
    xcrit=(((U/(K*(ka^m)))^(1/(h*m)))*((1/Sc)^(n/(h*m))));
    xcrit2=(((U/(K*Fms*(ka^m)))^(1/(h*m)))*((1/Sc)^(n/(h*m))));
    zcrit= (((U/K)^(1/n))*((ka^(-theta))/(1-(h*theta)))*(((L)^(1-(h*theta))) - (x_d((int16(((xcrit-xc)/dx))-1))^(1-(h*theta)))));
    zcrit2= (((U/(K*Fms))^(1/n))*((ka^(-theta))/(1-(h*theta)))*(((L)^(1-(h*theta))) - (x_d((int16(((xcrit-xc)/dx))-1))^(1-(h*theta)))));
    z0 = ((U/K)^(1/n))*((ka^-theta)/(1-h*theta)).*((x_d(end)^(1-h*theta)) - (x_d.^(1-h*theta)));
    dzdx0(1:Nnodes-1) = (1/dx)*(z0(1:Nnodes-1)-z0(2:Nnodes)); 
    z2 = ((U/(K*Fms))^(1/n))*((ka^-theta)/(1-h*theta)).*((x_d(end)^(1-h*theta)) - (x_d.^(1-h*theta)));
    dzdx2(1:Nnodes-1) = (1/dx)*(z2(1:Nnodes-1)-z2(2:Nnodes)); 
    if Fms<1
%         for i=1:int16((xcrit-xc+dx)/dx)
%             zwf(i)=z2(i)-z2(int16((xcrit-xc)/dx))+z0(int16((xcrit-xc)/dx));
%             dzdx_wf(i)=dzdx2(i);
%         end
%         for i =int16((xcrit-xc+dx)/dx):Nnodes
%             zwf(i)=z0(i);
%             dzdx_wf(i)=(dzdx0(i));
%         end
        load('D:\Users\srothman\Documents\fvd_tvd_slowwfssprofs.mat')
        load('D:\Users\srothman\Documents\ztvd_slowwfss.mat')
        zwf=z_slowwfss;
        dzdx_wf=dzdx_swfss;
    end
    if Fms>1
        xcrit_2=(((U/(K*Fms*(ka^m)))^(1/(h*m)))*((1/Sc)^(n/(h*m))));
        zcrit2=z0(int16((xcrit-xc)/dx))+(Sc*(xcrit-xcrit_2));
        for i=1:int16((xcrit_2-xc)/dx)
            zwf(i)=(z2(i)-z2(int16((xcrit_2-xc)/dx)) + z0(int16((xcrit-xc)/dx))+(Sc*(xcrit-xcrit_2)));
            dzdx_wf(i)=dzdx2(i);
        end
        for i =((int16((xcrit_2-xc)/dx)-1):int16((xcrit-xc)/dx))
            zwf(i)=(z0(int16((xcrit-xc)/dx))+((xcrit-x_d(i))*Sc));
            dzdx_wf(i)=Sc;
        end
        for i = int16(((xcrit-xc)/dx)): int16((L-xc)/dx)
            zwf(i)=z0(i);
            dzdx_wf(i)=dzdx0(i);
        end
        load('fvd_tvd_fastwfssprofs.mat')
        
    end
    
end

z=zwf;       %elevation for simple explicit upwind

%z=z0;
z(1:end-1) = z(1:end-1); %+Upulse;
zEX = z;    % elevation for Campforts solution
zTVD =z;

%% initial calculations

A = ka*x_d.^h; % calculate drainage area using Hack's law, which is formulated as a function of distance from the divide
A_sfo = flip(A);  %make first element of A be area at outlet
c = -K*A.^m;     % wave speed vector
c_ori=c;
c_m = min(0,c);
c_p = max(0,c);


%Campforts adds in ghost node at headwaters
% Extend rivers
xfvm=[x_d(1) x_d];

Afvm=[A(1) A];

% Extend inisurf

%zTVD = [zTVD(1) zTVD]; %UNUNCOMMENT IF not loading 
nb = length(xfvm);
s1 = [s1(1) s1];  %UNUNCOMMENT IF not loading 
if Fms<1
    s1=s1_slowwfss;
    zTVD=zpriorTVD_slowwfss;
end

%wave speed that starts from the divide
a = -K.*Afvm.^m;
a_ori=a;
a_m = min(0,a);
a_p = max(0,a);
a_vect=a_ori;
chi = dx * cumsum(flip(A).^(-theta)); %calculate chi at each stream node by integrating upstream from baselevel, which is why must flip A so A(1) is area at outlet

% check the Courant-Friedrichs-Lewy (CFL) Condition for stability
Amax = ka*L^h; %maximum drainage area (at outlet)  for mainstem
CFL = round((K*Amax^m)*(dt/dx),3);
disp(['CFL = ' num2str(CFL)])
if Fms>1
    load('fvd_tvd_fastwfssprofs.mat')
end
if (K*Amax^m)*(dt/dx)>1
    print('unstable parameters, need to increase dx or decrease dt')
end

%% run the model
for ti = 1:length(ta)
    
    if rem(ta(ti),plotevery)==0 %plot initial condition and then every plotevery iterations
        
        fh1 = figure(1);           %plot mainstem and tribs
        subplot(3,1,1)
        
        yyaxis left
         %plot(x_d/1e3,(z0),'r-')
         %hold on
        %plot(x_d/1e3,z, 'r-')
        %hold on
        ylim([0 (zmax+500)])
        xlim([0 L/1e3])
%         if ~isempty(wfix)
%             plot(x_d(wfix)/1e3,(z(wfix)),'r.')  
%         end
        
        %plot(x_d/1e3,zEX,'k-')
        %plot(x_d/1e3,(zTVD),'g-')
        
        if ~isempty(wfixTVD)
            plot(xfvm(wfixTVD)/1e3,(zTVD(wfixTVD)),'m.')  
        end
        hold on
        plot(xfvm/1e3,(zTVD),'g-')   %added  
        hold on
        plot(xcrit/1e3, zcrit, 'r*')
        hold on
        plot(xcrit2/1e3, zcrit2, 'g*')
        text(1,1000,['t = ' num2str(ta(ti)/1e6) ' Ma'])
        ylim([0 (zmax+500)])
        ylabel('Elevation (m)')
        xlabel('Distance from outlet (km)')
        title('FVM')
        hold off
        box on
        
        
        yyaxis right
        %plot(x_d/1e3,dzdx)
        plot(xfvm/1e3,s1)
        hold on
        plot(xcrit/1e3, Sc, 'r*')
        plot(xcrit2/1e3, Sc, 'g*')
        hold off
        ylim([0 1])
        
        
        subplot(3,1,2)
        
        yyaxis left
         %plot(x_d/1e3,(z0),'r-')
         %hold on
        plot(x_d/1e3,z, 'r-')
        hold on
        scatter(x_d(dzdx>Sc)/1e3,z(dzdx>Sc), 'm.')
        hold on
        plot(xcrit/1e3, zcrit, 'r*')
        plot(xcrit2/1e3, zcrit2, 'g*')

        ylim([0 (zmax+500)])
        xlim([0 L/1e3])
%         if ~isempty(wfix)
%             plot(x_d(wfix)/1e3,(z(wfix)),'r.')  
%         end
        title('FDM')
        
        %plot(x_d/1e3,zEX,'k-')
        %plot(x_d/1e3,(zTVD),'g-')
        
        %if ~isempty(wfixTVD)
        %    plot(x_d(wfixTVD)/1e3,(zTVD(wfixTVD)),'m.')  
        %end
        %hold on
        %plot(x_d/1e3,(zTVD),'g-')   %added     
        text(1,1000,['t = ' num2str(ta(ti)/1e6) ' Ma'])
        ylabel('Elevation (m)')
        xlabel('Distance from outlet (km)')
        hold off
        box on
        
        
        yyaxis right
        plot(x_d/1e3,dzdx)
        hold on
        plot(xcrit/1e3, Sc, 'r*')
        plot(xcrit2/1e3, Sc, 'g*')

        %plot(x_d/1e3,s1)
        hold off
        ylim([0 1])
        
        
        subplot(3,1,3)
        plot(flip(chi),z,'b')
        hold on
        set(gca, 'XDir','reverse')
        ylabel('Elevation (m)')
        xlabel('\chi (m)')
        xlim([0 max(chi)])
        ylim([0 zmax+500])
        hold off
        box on
        
%         subplot(3,1,3)
%         plot(x/1e3,dzdx,'k')
%         xlim([0 L/1e3])
%         %ylim([0 zmax])
%         ylabel('slope (m/m)')
%         xlabel('Distance from outlet (km)')
%         hold off
%         box on
        
        drawnow
    end
       
        
    if rem(ta(ti),saveevery)==0
        savedzdx_v(counter, : )=s1;
        savedzdx_d(counter, :)=dzdx;
        savez_v(counter, :)=zTVD;
        savez_d(counter, :)=z;
        savewf_v( counter, s1>Sc)=ones(1,sum(s1>Sc));
        savewf_d(counter, dzdx>Sc)=ones(1,sum(dzdx>Sc));
        %exportgraphics(fh1,['InstantaneouslongPulse_WFSteadyState_L100km_U200_WFT1F0p5_WFT2F2_' num2str(iteration) '.svg'],'ContentType','vector')
        %%SWMExportFig(gcf,['InstantaneousUniformPulse_UniformBackground_' num2str(iteration) '.eps'],'one-col',4.75,8,'-eps', '-RGB')
        counter = counter + 1;
        
       
    end
%% simple upwind FD
    %First uplift
    z(1:end-1) = z(1:end-1) + Ufield(1:end-1)*dt+ Pulseovertime(ti);
    
    %upwind difference, in this case with fix baselevel at x(end): right minus left differences because "upwind" side is the opposite side of the direction the wave is traveling, which is upstream.
    %note this never changes dzdx(0), which effectively fixes the slope at
    %baseleve =0 (the value we initialized it at above
    
    dzdx(1:Nnodes-1) = (1/dx)*(z(1:Nnodes-1)-z(2:Nnodes)); 
    bin_wfs_d=dzdx>Sc;
    counter_d=counter_d+bin_wfs_d;
    zold = z;
    E = K .* A.^m .* dzdx.^n;
    
    if waterfalls_ms==1%waterfalls_ms
        wfix = find(dzdx>Sc);                    %get index of nodes that have waterfalls
        E(dzdx>Sc) = Fms*E(dzdx>Sc);
    end
    
    z = z - (E*dt);   %update new elevations
    
%% Campforts explicit TVD-FVS  solution
    
     % %First uplift every point but boundary
    zTVD(1:end-1) = zTVD(1:end-1) + U*dt+ Pulseovertime(ti);   % z(1:end-1) = z(1:end-1) + Ufield(1:end-1)*dt+ Pulseovertime(ti);
    
    zpriorTVD=zTVD; 
    
    %calc slope even if n = 1
    TVD_r_t=[zTVD(2:end) zTVD(end)];
    s1=(zTVD-TVD_r_t)/dx;            %calc upwind slope
    
    %% TVD
    zTVD(zTVD<1e-6)=0;
    if n~=1
        TVD_r_t=[zTVD(2:end) zTVD(end)];
        s1=(zTVD-TVD_r_t)/dx;            %calc upwind slope
        s1(s1<1e-4)=0;
        exp_f=s1.^(n-1);
        exp_f(isinf(exp_f))=1;
        a=a_ori.*exp_f;
        a_m = min(0,a);
        a_p = max(0,a);
    end
    
    
    TVD_r=[zTVD(2:end) zTVD(end)];
    TVD_r=[zTVD(2:end) zTVD(end)];
    TVD_r2=[zTVD(3:end) zTVD(end) zTVD(end)];
    TVD_l=[zTVD(1) zTVD(1:end-1)];
    r_TVD=(TVD_r2-TVD_r)./(TVD_r-zTVD);
    r_TVD(diff(zTVD)==0)=1;
    r_TVD(1) = 1;
    %      r_TVD(nb-1) = r_TVD(nb-2);
    r_TVD(nb) = 1;
    
    %   Define Flux Limiter function
    %VANLEER
    phi = (r_TVD + abs(r_TVD))./(1 + abs(r_TVD));
    
    %   Compute fluxes for TVD
    TVD_r=[zTVD(2:end) zTVD(end)];
    TVD_l=[zTVD(1) zTVD(1:end-1)];
    
    F_rl = a_p.*zTVD + a_m.*TVD_r;
    F_rh = (1/2)*a.*(zTVD+TVD_r) - (1/2)*(a.^2).*(dt/dx).*(TVD_r-zTVD);
    F_ll = a_p.*TVD_l + a_m.*zTVD;
    F_lh= (1/2)*a.*(TVD_l+zTVD) - (1/2)*(a.^2).*(dt/dx).*(zTVD-TVD_l);
    %   Compute nz.ext time step
    phi_l=[phi(1) phi(1:end-1)];
    F_right = F_rl + phi.*(F_rh-F_rl);
    F_left = F_ll+ phi_l.*( F_lh- F_ll);
    ETVD = dt*(F_right-F_left)/dx;   %calc erosion
    
    if waterfalls_ms==1               %change erosion rate based on waterfall Fms
        wfixTVD = find(s1>Sc);        %get index of nodes that have waterfalls
        ETVD(s1>Sc) = Fms*ETVD(s1>Sc);
    end
    
    TVD_next= zTVD-ETVD;
    
    %   UPDATE info
    TVD_next(1) = TVD_next(2);
    TVD_next(nb) = zTVD(end);
    if any(~isreal(TVD_next))
        error('imaginary number');
        %         real(TVD_next);
    end
    zTVD = TVD_next;
    
    


  
end


%% save state 
%save('SteadyState_L100km_dx50_U50_Tilt.mat') %save the entire workspace to use
%as an IC

figure
   subplot(2,1,1)
   yyaxis left
   plot(x_d/1e3, dzdx)
   hold on
   plot(xcrit/1e3, Sc)
   plot(xcrit2/1e3, Sc)
   hold off
   ylim([0, .55])
   yyaxis right
   plot(x_d/1e3, counter_d)
   hold on
   xline(xcrit/1e3)
   xline(xcrit2/1e3)
   hold off
   ylim([-50000, 300000])
   xlim([0,20])

   subplot(2,1,2)
   yyaxis left
   plot(xfvm/1e3, s1)
   hold on
   plot(xcrit/1e3, Sc)
   plot(xcrit2/1e3, Sc)
   hold off
   ylim([0, .55])
   yyaxis right
   plot(x_d/1e3, counter_v)
   hold on
   xline(xcrit/1e3)
   xline(xcrit2/1e3)
   hold off
   ylim([-50000, 300000])
   xlim([0,20])

    %% 
   %load(fvm_fdm_fastwfss.mat)
   figure
   yyaxis left
   plot(x_d/1e3, dzdx, 'b', 'Linewidth', 5)
   hold on
   plot(xfvm/1e3, s1,  'g-')
   plot(xcrit/1e3, Sc)
   plot(xcrit2/1e3, Sc)
   hold off
   ylim([0, .55])
   ylabel('Slope')

   yyaxis right

   plot(x_d/1e3, counter_d, 'k', 'Linewidth', 5)
   
   hold on
   plot(x_d/1e3, counter_v, 'r-')
   %xline(xcrit/1e3)
   %xline(xcrit2/1e3)
   hold off
   ylabel('Waterfall Frequency')
   ylim([-10000, 120000])
   xlabel('Distance from Divide (km)')
   legend('Slope-FDM', 'Slope-FVM', 'xcritF1', 'xcritwf','WF Freq- FDM', 'WF Freq-FVM', 'Location', 'northeast')
   %xlim([0,20])

figure
sgtitle('Transient Response to 50 m Uplift over .5 MY')
subplot(2,1,1)
toplot=[5,14,24,35, 67];

for i=toplot
    semilogy(xfvm/1e3, savedzdx_v(i,:))
    hold on
    %legend('T = ' +num2str(i*saveevery) +' MY')
end
ylim([0, .55])
xlabel('Distance from Divide (km)')
ylabel('Slope')
title('Finite Volume Method')
ylim([1e-2, 4e-1])
hold off
subplot(2,1,2)

for i=toplot
    semilogy(x_d/1e3, savedzdx_d(i,:))
    hold on
end
xlabel('Distance from Divide (km)')
ylabel('Slope')
title("Finite Difference Method")
   ylim([0, .55])
legendyears=(toplot-1).*saveevery/1e6;
legend(string(legendyears(1)) +'Myr', string(legendyears(2)) +' Myr', string(legendyears(3))+' Myr', string(legendyears(4))+' Myr', string(legendyears(5)) +' Myr')
ylim([1e-2, 4e-1])
hold off

figure
sgtitle('Final Profiles after 50 m uplift pulse over .5 MY')
subplot(2,1,1)
plot(x_d/1e3, z_fastwfss, 'LineWidth', 5)
hold on
plot(xfvm/1e3, zpriorTVD_fastwfss, 'LineWidth', 2.5)
hold on
plot(x_d/1e3, z, 'LineWidth', 2.5)
hold on
plot(xfvm/1e3, zTVD, 'LineWidth', 1.5)
xlabel('Distance from Divide (km)')
ylabel('Elevation (m)')
hold off
legend('pre-uplift', 'post uplift FDM', 'post uplift FVM')

subplot(2,1,2)
semilogy(x_d/1e3, dzdx_fastwfss,'LineWidth', 5)
hold on
semilogy(xfvm/1e3, s1_fastwfss,'LineWidth', 1.5)
hold on

semilogy(x_d/1e3, dzdx,'LineWidth', 2.5)
hold on
semilogy(xfvm/1e3, s1,'LineWidth', 1.5)
   ylim([0, .55])
xlabel('Distance from Divide (km)')
ylabel('Slope')
hold off


%% j
figure
semilogy(xfvm/1e3, s1_fastwfss,'LineWidth', 1.5)


s1_fastwfss(1:Nnodes) = (1/dx)*(zpriorTVD_fastwfss(1:Nnodes)-zpriorTVD_fastwfss(2:Nnodes+1)); 