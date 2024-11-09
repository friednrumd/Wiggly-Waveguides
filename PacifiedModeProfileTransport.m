%% Initialize

clc
close all
clearvars -except rawmodes rawfoms xy dA psi0 N nmodes n0 eps

addpath('C:\Users\natef\OneDrive - University of Maryland\MATLAB\DataFiles')
addpath('.\Functions')


c_const=299792458;
mu_const=4*pi*10^(-7);
eps_const=1/(c_const^2*mu_const);
ell=1.55/2/pi;                      %length unit

if exist('rawmodes','var') == 0
    rawmodes    = readmatrix("ModeProfile3DProfileScaled.txt");

    rawfoms     = readmatrix("ModeProfileFOMRestricted.txt");

    xy=rawmodes(:,3:4);
    dA=rawmodes(:,6);
    eps=rawmodes(:,7);

    N       =size(rawmodes  ,1);
    nmodes  =size(rawfoms,1);
    
    INVERT=0;

    psi0=zeros(N,3,2,nmodes);

    for i=1:nmodes
        for j=1:2
            psi0(:,j,1,i)=  rawmodes(:,2+11*(i-1)+j-1+6)/sqrt((rawfoms(i,end)));
            psi0(:,j,2,i)=  (-1)^INVERT*rawmodes(:,2+11*(i-1)+j-1+9)/sqrt((rawfoms(i,end)));
        end
            psi0(:,3,1,i)=  (-1)^INVERT*rawmodes(:,2+11*(i-1)+3-1+6)/sqrt((rawfoms(i,end)));
            psi0(:,3,2,i)=  rawmodes(:,2+11*(i-1)+3-1+9)/sqrt((rawfoms(i,end)));
    end

    psi0=flip(psi0,4);
    rawfoms=flip(rawfoms,1);
    n0      =rawfoms(:,3);
end

FOMs=["nc", "Complex Index", "Real Index", "Imag Index", "Core Power", "Clad Power", "PML Power", "Sub Power", "Air Power", "Total Power", "Core S", "Clad S", "PML S", "Sub g", "Air S", "Total S", "Core power fraction", "Ey/Ex (rad)", "TE/TM (rad)", "Total z Power"];


%% Settings
dwTotal=0.100       /ell;  %um
divs=(1:0.25:4);        %orders of magnitudes of fractions of total change to compute by (dw=dwTotal./10.^divs)

MODES=1:8;
nmodes=length(MODES);


%% Parameters
ecore=2.1482;
eclad=3.9851;


corewidth=3.000     /ell;
coreheight=0.400    /ell;

cladwidth=18.000     /ell;
cladheight=(3.55+2)     /ell;

pmlt=2*1.55 /ell;

Atotal=(2*pmlt+cladheight)*(2*pmlt+cladwidth);

%% Preparations

w0=corewidth/2;
dw=dwTotal./10.^(divs);

psi0=psi0(:,:,:,MODES);

%% Normalization

psi=zeros(N,3,2,nmodes);
for i=1:nmodes
    psi(:,:,:,i)=psi0(:,:,:,i)/sqrt(conjS_Metric(psi0(:,:,:,i),psi0(:,:,:,i),dA));
end


%% Integrals

S0=zeros(nmodes,nmodes);
H0=zeros(nmodes,nmodes);
T0=zeros(nmodes,nmodes);

for i=1:nmodes
    S0(i,i)=S_Metric(psi(:,:,:,i),psi(:,:,:,i),dA);
    H0(i,i)=H_Operator(psi(:,:,:,i),psi(:,:,:,i),dA,eps);
    T0(i,i)=1/2*( ( psi(:,3,1,i) )'*(dA.*eps.*psi(:,3,1,i))+( psi(:,3,2,i) )'*( dA.*psi(:,3,2,i)) );

    for j=(1+i):nmodes
        S0(i,j)=S_Metric(psi(:,:,:,i),psi(:,:,:,j),dA);
        H0(i,j)=H_Operator(psi(:,:,:,i),psi(:,:,:,j),dA,eps);
        T0(i,j)=1/2*( ( psi(:,3,1,i) )'*(dA.*eps.*psi(:,3,1,j))+( psi(:,3,2,i) )'*( dA.*psi(:,3,2,j)) );

        S0(j,i)=(S0(i,j));
        H0(j,i)=(H0(i,j));
        T0(j,i)=(T0(i,j));

    end
end

V0=H0-T0;


%% Coefficient Matrices
cmap=readmatrix('RedBlueColormap.txt');

figure
imagesc(MODES,MODES,real(S0))
colorbar
colormap(cmap)
axis xy
title('S0')
clim([-2,2])

figure
imagesc(MODES,MODES,real(H0))
colorbar
colormap(cmap)
axis xy
title('H0')
clim([-2,2])


figure
imagesc(MODES,MODES,real(T0))
colorbar
colormap(cmap)
axis xy
title('T0')
clim([-2,2])


figure
imagesc(MODES,MODES,real(V0))
colorbar
colormap(cmap)
axis xy
title('V0')
clim([-2,2])






%% Running

% for d=divs
%     modalbasisprofiles=zeros(nmodes,nmodes,d);
%     eigs=zeros(nmodes,d);
% 
%     modalbasisprofiles(:,:,1)=eye(nmodes);
%     eigs(:,1)=n0;
% 
%     cw=0;
%     for w=w0:dWtotal/d:(w0+dwTotal)
%         cw=cw+1;
% 
% 
%     end
% 
% 

%end
%% Post Anaylysis


%% DIAGNOSTICS



%% Mode Fields

mode1=1;
mode2=2;
area=10;

max1=(max(max(max(abs(psi0(:,:,:,mode1))))));
max2=(max(max(max(abs(psi0(:,:,:,mode2))))));

fields=["E_x"; "E_y"; "E_z"; "B_x"; "B_y"; "B_z"];
figure
for m=1:3
    subplot(2,3,m)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*real(psi0(:,m,1,mode1)))
    colorbar
    title(fields(m))
    xlabel('x (\mum)')
    ylabel('y (\mum)')

    subplot(2,3,m+3)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*real(psi0(:,m,2,mode1)))
    colorbar
    title(fields(m+3))
    xlabel('x (\mum)')
    ylabel('y (\mum)')

end
sgtitle('Real Mode Fields')

figure
for m=1:3
    subplot(2,3,m)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*imag(psi0(:,m,1,mode1)))
    colorbar
    title(fields(m))
    xlabel('x (\mum)')
    ylabel('y (\mum)')


    subplot(2,3,m+3)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*imag(psi0(:,m,2,mode1)))

    colorbar
    title(fields(m+3))
    xlabel('x (\mum)')
    ylabel('y (\mum)')

end
sgtitle('Imag Mode Fields')

%% Particular Mode Field

figure
scatter(xy(:,1),xy(:,2),area*dA,abs(psi0(:,3,2,5)))
colorbar
title(['Mode' num2str(mode1)])
xlabel('x (\mum)')
ylabel('y (\mum)')

%% FOM Plots

for i=3:(size(rawfoms,2))

    figure
    plot((MODES),real(rawfoms(MODES,i)),'.r','MarkerSize',10)
    title(FOMs(i))
    xlabel('Mode #')

end





% 
%  %% colormap
% radius=50;
% 
% cmap=zeros(3,2*radius+1);
% cmap(:,radius+1)=[1; 1; 1];
% 
% for i=1:radius
%     cmap(:,radius+1-i)=[1 1-(i)/radius 1-(i)/radius ];
%     cmap(:,radius+1+i)=[1-(i)/radius, 1-(i)/radius, 1];
% end
% 
% fprintf(fopen('RedBlueColormap.txt','w'),'%i %i %i \n',flip(cmap));
% 
% %% read
% clc
% clear variables
% 
% 
% q=readmatrix('RedBlueColormap.txt');







