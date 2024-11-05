%% Initialize

clc
close all
clearvars -except modedata

addpath('.\Functions')
addpath('.\DataFiles')

if exist('modedata') == 0
    modedata = table2array(readtable("C:\Users\natef\OneDrive - University of Maryland\MATLAB\Floquet Research\DataFiles\poorcrosssectionmodes3d.txt"));

    i=1;
    while modedata(i,1)==modedata(1,1)
        i=i+1;
    end
    N(1)=i-1;
    
    if rem(size(modedata,1),N(1))==0
        N(2)=size(modedata,1)/N(1);
    else
        fprintf('Ny \n');
        scream()
        return
    end

    nmodes=round((size(modedata(1,:),2)-2)/3,1);
    N(3)=N(1)*N(2);
    
    ModeData=zeros(nmodes,N(1),N(2));

    for j=1:nmodes
        ModeData(j,1:N(1),1:N(2))=reshape(modedata(:,2+j),[N(1),N(2)]);
    end

    Ex0=zeros(N(1),N(2),nmodes);
    Ey0=zeros(N(1),N(2),nmodes);
    Ez0=zeros(N(1),N(2),nmodes);
    for i=1:nmodes
        Ex0(:,nmodes)=modedata(:,3*i);
        Ey0(:,nmodes)=modedata(:,3*i+1);
        Ez0(:,nmodes)=modedata(:,3*i+2);
    end
end

n0=[3.437208 3.408583 3.360361 3.291731 3.201465 3.087806 2.948290 2.779468 2.576498 2.332719 2.040777 1.713245];


%% Settings
ell=1.55/2/pi; %length unit

dwTotal=0.100    /ell; %um

divs=(1:0.25:4); %fractions of total change to compute by

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

dA=(xVals(N(1)+1,1)-xVals(1,1))*(xVals(2,2)-xVals(1,2));

coref=(heaviside(xVals(:,1)'+corewidth/2)-heaviside(xVals(:,1)'-corewidth/2)).*(heaviside(xVals(:,2)+cladwidth/2)-heaviside(xVals(:,2)-cladwidth/2));

H0=zeros(nmodes,nmodes);
for i=1:nmodes
    for j=1:nmodes
        H0=( (Ex0(:,:))'*coref*Ex0(:,:)+(Ey0(:,:))'*coref*Ey0(:,:)+(Ez0(:,:))'*coref*Ez0(:,:) )*dA;
    end
end

S0=zeros(nmodes,nmodes);


imagesc(1:nmodes,1:nmodes,log(abs(H0)))
colorbar
colormap('jet')
axis xy


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


