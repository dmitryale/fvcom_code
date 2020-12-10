% function[]=get_siglocs();
% call from Massimo3.m; L196
% Testing interp_azonal3D_v2 06/01/2017
% (c) dmitry.aleynik@sams.ac.uk 2017.07.21;                               &
%  __o_O__¬                                                               &
%  \_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS                               &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 global Mesh %CVM;
% load('Mesh.mat');
% load('CVM_090101.mat');

a1u=Mesh.a1u;
a2u=Mesh.a2u;
nbe=Mesh.nbe;
  xc=Mesh.uvnode(:,1);
  yc=Mesh.uvnode(:,2);
   x=Mesh.nodexy(:,1);
   y=Mesh.nodexy(:,2);

% Bathymetry
n=Mesh.Nelems; %n=79244;
%  bath=100*(-abs(randn(n,1)));

% Sigma layers
% sigloc=[0,-0.0219698204412441,-0.0676189020719762,-0.155592441548392,...
%     -0.302936674163750,-0.5,-0.697063325836250,-0.844407558451608,...
%     -0.932381097928024,-0.978030179558756,-1]';

%%!da---- see 'sigma.dat' : SIGMA COORDINATE TYPE = TANH;  NUMBER OF SIGMA LEVELS = 11
    kbm=11; kbm1=kbm-1;
    kb =11; %!da correct
       DU2 = 2.0;  DL2 = 2.0; clear siglay siglev
 for K=1:kbm1,...
     X1=DL2+DU2;
     X1=X1*(kbm1-K)/kbm1;
     X1=X1-DL2;
     X1=tanh(X1);
     X2=tanh(DL2);
     X3=X2+tanh(DU2);
   siglev(K+1,1)=(X1+X2)/X3 - 1.0;
 end;
 for K=1:kbm1,...
    siglay(K,1) = 0.5*(siglev(K)+siglev(K+1));
 end
%     siglay(kb,1)  = 2.0*siglay(kb-1)-siglay(kb-2); % no need 11th
%% !da------------------------
   sigloc=siglev;
% Triangles all
tri=Mesh.trinodes;
    Mesh.sigloc=sigloc;
    Mesh.siglev=siglev;%11
    Mesh.siglay=siglay;%10
 fileToRead1='sigma.dat';
      zl=[0. 1.0  5.0 10.0 20.0 30.0      50.0     100.      150. 200.  225.];  % 11
   [ z,z1,zz,zz1,dz,dz1,dzz1,SIGMA] = setup_sigma_gen(fileToRead1,zl,Mesh );
    Mesh.z1  = z1;
    Mesh.dz1 = dz1;
    Mesh.zz1 = zz1;
    Mesh.dzz1= dzz1;


%%
