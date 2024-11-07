%% Initialize

clc
close all
clearvars -except rawmodes rawindices xy dA psi0 N nmodes

addpath('C:\Users\natef\OneDrive - University of Maryland\MATLAB\DataFiles')
addpath('.\Functions')

c_const=299792458;
mu_const=4*pi*10^(-7);
eps_const=1/(c_const^2*mu_const);
ell=1.55/2/pi;                      %length unit

if exist('rawmodes','var') == 0
    rawmodes = readmatrix("ModeProfile3DProfile.txt")/eps_const/10^12;
    rawindices = readmatrix("ModeProfileFOMRestricted.txt");
    xy=rawmodes(:,1:2)/ell;
    dA=rawmodes(:,3)*10^12/ell^2;

    N       =size(rawmodes  ,1);
    nmodes  =size(rawindices,1);

    psi0=zeros(N,3,2,nmodes);
    for i=1:nmodes
        for j=1:3
            psi0(:,j,1,i)=  rawmodes(:,1+6*(i-1)+(j-1)+3);
            psi0(:,j,2,i)=  rawmodes(:,1+6*(i-1)+(j-1)+6)*c_const;
        end
    end
end


%% Settings
dwTotal=0.100       /ell;  %um
divs=(1:0.25:4);        %orders of magnitudes of fractions of total change to compute by (dw=dwTotal./10.^divs)

TRUNCATE_MODES=1:8;



%% Parameters
ncore=2.03-0.1*1i;
nclad=1.46-0.05*1i;

corewidth=3.000     /ell;
cladwidth=9.000     /ell;

coreheight=0.400    /ell;
cladheight=3.000    /ell;

%% Preparations

ecore=ncore^2;
eclad=nclad^2;

w0=corewidth/2;
dw=dwTotal./10.^(divs);

coref= @(r) (ecore-eclad)*(heaviside(r(:,1)+w0)-heaviside(r(:,1)-w0)).*(heaviside(r(:,2)+coreheight/2)-heaviside(r(:,2)-coreheight/2))+eclad;

nmodes=length(TRUNCATE_MODES);
psi0=psi0(:,:,:,TRUNCATE_MODES);

S0=zeros(nmodes,nmodes);
H0=zeros(nmodes,nmodes);
for i=1:nmodes
    S0(i,i)=    S_Metric(psi0(:,:,:,i),psi0(:,:,:,i),dA);
    H0(i,i)=    S_Metric(psi0(:,:,:,i),psi0(:,:,:,i),dA);

    for j=i+1:nmodes
        S0(i,j)=S_Metric(psi0(:,:,:,i),psi0(:,:,:,j),dA);
        S0(j,i)=S0(i,j);

        H0(i,j)=S_Metric(psi0(:,:,:,i),psi0(:,:,:,j),dA);
        H0(j,i)=H0(i,j);
    end
end
modes=1:8;

figure
imagesc(modes,modes,(abs(S0(modes,modes))))
colorbar
colormap('jet')
axis xy
title('S0 Arg')

figure
imagesc(modes,modes,(angle(S0(modes,modes))))
colorbar
colormap('jet')
axis xy
title('S0 Angle')

q1=abs(S0);
% Normalize
for i=1:nmodes
    S0(i,i)=    1;
    H0(i,i)=    H0(i,i)/S0(i,i);
    psi0(:,:,:,i)=psi0(:,:,:,i)/sqrt(S0(i,i));

    for j=i+1:nmodes
        S0(i,j)=S0(i,j)/sqrt(S0(i,i)*S0(j,j));
        S0(j,i)=S0(i,j);

        H0(i,j)=H0(i,j)/sqrt(S0(i,i)*S0(j,j));
        H0(j,i)=H0(i,j);
    end
end

figure
imagesc(modes,modes,(abs(S0(modes,modes))))
colorbar
colormap('jet')
axis xy
title('Normalized S0 Arg')

figure
imagesc(modes,modes,(angle(S0(modes,modes))))
colorbar
colormap('jet')
axis xy
title('Normalized S0 Angle')
q2=abs(S0);


%% Running

for d=divs
    modalbasisprofiles=zeros(nmodes,nmodes,d);
    eigs=zeros(nmodes,d);

    modalbasisprofiles(:,:,1)=eye(nmodes);
    eigs(:,1)=n0;
    
    cw=0;
    for w=w0:dWtotal/d:(w0+dwTotal)
        cw=cw+1;


    end



end
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
    clim([-max1,max1])

    subplot(2,3,m+3)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*real(psi0(:,m,2,mode1)))
    colorbar
    title(fields(m+3))
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    clim([-max1,max1])

end
sgtitle('Real Mode Fields')

figure
for m=1:3
    subplot(2,3,m)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*real(psi0(:,m,1,mode1)))
    colorbar
    title(fields(m))
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    clim([-max1,max1])


    subplot(2,3,m+3)
    scatter(xy(:,1),xy(:,2),area*dA,area*dA.*real(psi0(:,m,2,mode1)))

    colorbar
    title(fields(m+3))
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    clim([-max1,max1])

end
sgtitle('Imag Mode Fields')

%% Particular Mode Field

figure
scatter(xy(:,1),xy(:,2),area*dA,area*dA.*imag(psi0(:,2,1,mode1)))
colorbar
title(['Mode' num2str(mode1)])
xlabel('x (\mum)')
ylabel('y (\mum)')

%% 



