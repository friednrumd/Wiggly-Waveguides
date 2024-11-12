%% Initialize

clc
close all
clearvars -except rawmodes rawfoms xy dA Psi_Init N nmodes n0 eps

addpath('C:\Users\natef\OneDrive - University of Maryland\MATLAB\DataFiles')
addpath('.\Functions')


c_const=299792458;
mu_const=4*pi*10^(-7);
eps_const=1/(c_const^2*mu_const);
ell=1.55/2/pi;                      %length unit

if exist('rawmodes','var') == 0

    rawmodes    = readmatrix("ModeProfile2DdataNorm.txt");
    rawfoms     = readmatrix("ModeProfileFOMs.txt");
    rawintegrals     = readmatrix("ModeProfileFOMs.txt");

    xy=rawmodes(:,3:4);
    dA=rawmodes(:,6);
    eps=rawmodes(:,7);

    N       =size(rawmodes  ,1);
    nmodes  =size(rawfoms,1);
    
    INVERT=0;

    Psi_Init=zeros(N,3,2,nmodes);

    for i=1:nmodes
        for j=1:2
            Psi_Init(:,j,1,i)=              rawmodes(:,2+12*(i-1)+(j-1)+7)*1i;
            Psi_Init(:,j,2,i)=  (-1)^INVERT*rawmodes(:,2+12*(i-1)+(j-1)+10)*1i;
        end
            Psi_Init(:,3,1,i)=  (-1)^INVERT*rawmodes(:,2+12*(i-1)+(3-1)+7)*1i;
            Psi_Init(:,3,2,i)=              rawmodes(:,2+12*(i-1)+(3-1)+10)*1i;
    end

    Psi_Init=flip(Psi_Init,4);
    rawfoms=flip(rawfoms,1);
    n0      =rawfoms(:,1);
end

FOMs=["nc", "Complex Index", "Real Index", "Imag Index", "Core Power", "Clad Power", "PML Power", "Sub Power", "Air Power", "Total Power", "Core S", "Clad S", "PML S", "Sub g", "Air S", "Total S", "Core power fraction", "Ey/Ex (rad)", "TE/TM (rad)", "Total z Power"];


%% Settings
dwTotal=0.05       /ell;  %um
divs=(0:0.20:6);        %orders of magnitudes of fractions of total change to compute by (dw=dwTotal./10.^divs)
divs=round(10.^(divs));

MODES=1:nmodes;


%% Parameters
ecore=2.1482;
eclad=3.9851;

corewidth=3.000     /ell;
coreheight=0.400    /ell;

cladwidth=18.000     /ell;
cladheight=(3.55+2)     /ell;

%% Preparations

nmodes=length(MODES);

w0=corewidth/2;
dw=dwTotal./(divs);

%% Normalization

Psi_0=zeros(N,3,2,nmodes);
for i=1:nmodes
    Psi_0(:,:,:,i)=Psi_Init(:,:,:,i)/sqrt((S_Metric((Psi_Init(:,:,:,i)),Psi_Init(:,:,:,i),dA)));
end

%% Integrals

S0 =zeros(nmodes,nmodes);
Hx0=zeros(nmodes,nmodes);
Hy0=zeros(nmodes,nmodes);
Hz0=zeros(nmodes,nmodes);


for i=1:nmodes

     S0(i,i)=  S_Metric((Psi_0(:,:,:,i)),Psi_0(:,:,:,i),dA);
    Hx0(i,i)=H_Operator((Psi_0(:,:,:,i)),Psi_0(:,:,:,i),dA,eps,1);
    Hy0(i,i)=H_Operator((Psi_0(:,:,:,i)),Psi_0(:,:,:,i),dA,eps,2);
    Hz0(i,i)=H_Operator((Psi_0(:,:,:,i)),Psi_0(:,:,:,i),dA,eps,3);

    for j=(1+i):nmodes

         S0(i,j)=  S_Metric((Psi_0(:,:,:,i)),Psi_0(:,:,:,j),dA);
        Hx0(i,j)=H_Operator((Psi_0(:,:,:,i)),Psi_0(:,:,:,j),dA,eps,1);
        Hy0(i,j)=H_Operator((Psi_0(:,:,:,i)),Psi_0(:,:,:,j),dA,eps,2);
        Hz0(i,j)=H_Operator((Psi_0(:,:,:,i)),Psi_0(:,:,:,j),dA,eps,3);

         S0(j,i)=( S0(i,j));
        Hx0(j,i)=(Hx0(i,j));
        Hy0(j,i)=(Hy0(i,j));
        Hz0(j,i)=(Hz0(i,j));
        
    end
end

H0=Hx0+Hy0+Hz0;
V0=Hx0+Hy0;

%% Coefficient Matrices
cmap=readmatrix('Murica.txt');

figure
imagesc(MODES,MODES,real(S0))
colorbar
colormap(cmap)
axis xy
title('S0')
clim([-2,2])


% figure
% imagesc(MODES,MODES,abs(Hx0))
% colorbar
% colormap(cmap)
% axis xy
% title('Hx0')
% clim([-2,2])
% 
% figure
% imagesc(MODES,MODES,abs(Hy0))
% colorbar
% colormap(cmap)
% axis xy
% title('Hy0')
% clim([-2,2])
% 
% figure
% imagesc(MODES,MODES,abs(Hz0))
% colorbar
% colormap(cmap)
% axis xy
% title('Hz0')
% clim([-2,2])

figure
imagesc(MODES,MODES,real(H0))
colorbar
colormap(cmap)
axis xy
title('H0')
clim([-2,2])


figure
plot(real(n0(MODES)),-imag(n0(MODES)),'r.',real(diag(H0)./diag(S0)),-imag(diag(H0)./diag(S0)),'.b','MarkerSize',20)
legend('From COMSOL','From Integrals','Location','NorthEast')
title('n-ik')
xlabel('n')
ylabel('k')






%% Transport

nurFinal=zeros(nmodes,nmodes,length(divs));
cd=0;
for dd=divs
    cd=cd+1;

    
    % figure
    % hold on
    % colorbar
    % colormap(cmap)
    % clim([-30,30])
    % 
    nur   =zeros(nmodes,nmodes,dd);
    nr     =zeros(nmodes,dd);
    
    nur(:,:,1)=eye(nmodes);
    nr(:,1)    =n0;
    
    cw=0;
    for w=ws
        cw=cw+1;
    
        % if mod(cw,100)==0
        %     imagesc(1:nmodes,1:nmodes,10*log10(abs(nur(:,:,cw))))
        %     title(['w= ' num2str(ws(cw)*ell)])
        %     drawnow
        %     clc
        %     fprintf('w=%1.3f which is %i percent done\n',ws(cw)*ell, ceil(100*cw/length(ws)))
        % end
    
    
        Hr=(nur(:,:,cw).')*(-w^(-2)*(Hz0+Hx0)+Hy0)*(nur(:,:,cw));
        gammar=(1-eye(nmodes))./(nr(:,cw)-nr(:,cw)'+10^-20).*Hr;
    
        nr(:,cw+1)=nr(:,cw)+diag(Hr)*dw(cd);
        nur(:,:,cw+1)=(eye(nmodes)-gammar*dw(cd))*nur(:,:,cw);
    
        % for i=1:nmodes
        %     nr(i,cw+1)=nr(i,cw)+(nur(i,:,cw))*(-w^(-2)*(Hz0+Hx0)+Hy0)*(nur(:,i,cw))*dwTotal/dd;
        %     nur(i,i,cw+1)=nur(i,i,cw);
        % 
        %     for j=i+1:nmodes
        %         nur(i,j,cw+1)=nur(i,j,cw)-1/(nr(i,cw)-nr(j,cw))*(nur(i,:,cw))*(-w^(-2)*(Hz0+Hx0)+Hy0)*(nur(:,j,cw))*dwTotal/dd;
        %         nur(j,i,cw+1)=-nur(i,j,cw);
        %     end
        % end
    
    end
    clc
    fprintf('dw=%0.4f nm, Delta w=%1.2f um, N= %i, at step %i of %i \n',dw(cd)*ell*1000,dwTotal*ell,dd, cd, length(divs))
    ws=w0:dw(cd):(w0+dwTotal);

    ws=[ws w+dwTotal/dd];
    
    nurFinal(:,:,cd)=nur(:,:,cw+1);
end

%% Sumsqr mode profile error


errornur=zeros(length(divs)-1,1);
for k=1:(length(divs)-1)
    q=(nurFinal(:,:,k)-nurFinal(:,:,end));
    errornur(k)=trace(q'*q);
end

loglog(dw(1:(end-1)),4700*(dw(1:(end-1))).^(2.09),'r','LineWidth',2)
hold on
loglog(dw(1:(end-1)),errornur,'.b','MarkerSize',20)
xlabel('dw (\mu m)')
ylabel('MSE compared with final profile')
ylim([0,1])

%end
%% Indices Plot

shortl=cw;
figure
for i=1:nmodes
    semilogy(ws(1:shortl)*ell-w0*ell,abs(nr(i,1:shortl)-nr(end,1:shortl)),'LineWidth',2)
    if i==1
    hold on
    end
end
xlabel('\Deltaw (\mum)')
ylabel('n-n_{min}')

%% nu animation

figure
hold on
for c=1:50:cw
    imagesc(1:nmodes,1:nmodes,10*log10(abs(nur(:,:,c))))
    colorbar
    colormap(cmap)
    clim([-30,30])
    title(['w= ' num2str(ws(c)*ell)])
    drawnow
    clc
    fprintf('w=%1.3f which is %i percent done\n',ws(c)*ell, round(100*c/cw))
    
end



%% Mode Fields

mode1=1;
area=10;


fields=["E_x"; "E_y"; "E_z"; "B_x"; "B_y"; "B_z"];
figure
for m=1:2
    subplot(2,3,m)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*real(Psi_0(:,m,1,mode1)))
    colorbar
    title(fields(m))
    xlabel('x (\mum)')
    ylabel('y (\mum)')

    subplot(2,3,m+3)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*real(Psi_0(:,m,2,mode1)))
    colorbar
    title(fields(m+3))
    xlabel('x (\mum)')
    ylabel('y (\mum)')

end
subplot(2,3,3)
scatter(xy(:,1),xy(:,2),area*dA,area*dA.*imag(Psi_0(:,3,1,mode1)))
colorbar
title(fields(3))
xlabel('x (\mum)')
ylabel('y (\mum)')

subplot(2,3,6)
scatter(xy(:,1),xy(:,2),area*dA,area*dA.*imag(Psi_0(:,3,2,mode1)))
colorbar
title(fields(6))
xlabel('x (\mum)')
ylabel('y (\mum)')
sgtitle('Primary Mode Field Components')

figure
for m=1:3
    subplot(2,3,m)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*angle(Psi_0(:,m,1,mode1)))
    colorbar
    title(fields(m))
    xlabel('x (\mum)')
    ylabel('y (\mum)')


    subplot(2,3,m+3)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*angle(Psi_0(:,m,2,mode1)))

    colorbar
    title(fields(m+3))
    xlabel('x (\mum)')
    ylabel('y (\mum)')

end
sgtitle('Mode Field Angles')

%% Particular Mode Field Plot
mode1=1;

figure
scatter(xy(:,1),xy(:,2),area*dA,real(Psi_0(:,1,1,mode1).*Psi_0(:,2,2,mode1)-Psi_0(:,1,2,mode1).*Psi_0(:,2,1,mode1)))
colorbar
title(['Mode' num2str(mode1) 'Power'])
xlabel('x (\mum)')
ylabel('y (\mum)')

%% FOM Plots

for i=3:(size(rawfoms,2))

    figure
    plot((MODES),real(rawfoms(MODES,i)),'.r','MarkerSize',10)
    title(FOMs(i))
    xlabel('Mode #')

end







