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

divs=flip(round(10.^(   (3:0.1:5)   )));       %orders of magnitudes of fractions of total change to compute by
MODES=1:nmodes;


%% Parameters
ecore=2.1482;
eclad=3.9851;

corewidth=3.000     /ell;
coreheight=0.400    /ell;

cladwidth=18.000     /ell;
cladheight=(3.55+2)     /ell;


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

w0=corewidth/2;
dw=dwTotal./(divs);

nurFinal=zeros(nmodes,nmodes,length(divs));
errornur=zeros(length(divs),1);

cd=0;
for dd=divs
    cd=cd+1;

    ws=w0:dw(cd):(w0+dwTotal);

    

    % 
    nur   =zeros(nmodes,nmodes,dd);
    nr     =zeros(nmodes,dd);
    
    nur(:,:,1)=eye(nmodes);
    nr(:,1)    =n0;
    
    cw=0;
    for w=ws
        cw=cw+1;
    
        % Hr=(nur(:,:,cw).')*(-(Hz0+Hx0)*w^(-2)+Hy0)*(nur(:,:,cw));
        % gammar=((1-eye(nmodes))./(nr(:,cw)-nr(:,cw)'+10^-20)).*Hr;
        % 
        % nr(:,cw+1)=nr(:,cw)+diag(Hr)*dw(cd);
        % nur(:,:,cw+1)=(eye(nmodes)-gammar*dw(cd))*nur(:,:,cw);
        
        % [dn1, dnu1]  =ParameterUpdateFirst(nr(:,cw), nur(:,:,cw),w,Hz0+Hx0, Hy0);
        % nr(:,cw+1)   =nr(:,cw)+dn1*dw(cd);
        % nur(:,:,cw+1)=nur(:,:,cw)+dnu1*dw(cd);
        
        [nrk1, nurk1]=ParameterUpdateFirst(nr(:,cw)              , nur(:,:,cw)               , w         ,Hz0+Hx0, Hy0);
        [nrk2, nurk2]=ParameterUpdateFirst(nr(:,cw)+nrk1*dw(cd)/2, nur(:,:,cw)+nurk1*dw(cd)/2, w+dw(cd)/2,Hz0+Hx0, Hy0);
        [nrk3, nurk3]=ParameterUpdateFirst(nr(:,cw)+nrk2*dw(cd)/2, nur(:,:,cw)+nurk2*dw(cd)/2, w+dw(cd)/2,Hz0+Hx0, Hy0);
        [nrk4, nurk4]=ParameterUpdateFirst(nr(:,cw)+nrk3*dw(cd)  , nur(:,:,cw)+nurk3*dw(cd)  , w+dw(cd)  ,Hz0+Hx0, Hy0);

        nr(:,cw+1)=   nr(:,cw)+dw(cd)/6*( nrk1+ 2*nrk2+ 2*nrk3+ nrk4);
        nur(:,:,cw+1)=nur(:,:,cw)+dw(cd)/6*(nurk1 +2*nurk2+2*nurk3+nurk4);
    end
    nurFinal(:,:,cd)=nur(:,:,cw+1);
    errornur(cd)=sqrt(sumsqr(abs(nurFinal(:,:,cd)-nurFinal(:,:,1))));

    clc
    fprintf('To transport by Delta w=%1.2f um, via %i steps of %f nm, the RMSE is about %1.2e. \nThe calculations are approximately %2.0f percent done \n',dwTotal*ell,dd,dw(cd)*ell*1000, errornur(cd), sum(divs(1:cd))/sum(divs)*100)
    %ws=[ws w+dwTotal/dd];
    
end

%% Sumsqr mode profile error

loglog(dw(1:(cd)),sqrt(4700*(dw(1:(cd))).^(2.09)),'r','LineWidth',2)
hold on
loglog(dw(1:(cd)),sqrt(errornur),'.b','MarkerSize',20)
xlabel('dw (\mu m)')
ylabel('RMSE')
legend('First Order Fit for \Deltaw=0.05 \mum','rk4')

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







