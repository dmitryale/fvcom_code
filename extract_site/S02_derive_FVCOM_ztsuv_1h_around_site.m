% function[] = S02_derive_FVCOM_ztsuv_1h_around_site ();
% use a cluster to get annual subsets (BC1): 
% > sbatch generic_job_Scope_matlab.sh
% > matlab -softwareopengl -r "S02_derive_FVCOM_ztsuv_1h_around_site,exit"&
% ========================================================================&
% To derive tidal (1h) FVCOM model currents data arround target location(s)
% INput: mtb,mte L:74; L:91;['/mi*nestnodes.dat'];['minch2_' 'Scope' '_nest.mat'];
% OUTput:[path_out 'BC1' num2str(ki) '_' cdates '.mat'];% L:357           &
%                                                                         &
% History : based on R05u_mk_residual_meanflow_1h_NL.m ; 2017.07.21       &
% (c) created : Dmitry Aleynik AT sams.ac.uk, 29/10/2018                  & 
%  __o_O__/X/                                                             &
%  \_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS, COMPASS/CAMPUS 2018-2021     &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

clear all
close all
       this_dir=pwd;  DRV=this_dir(1:1);
 mk_newf=0; % to send into mk_figures_integral.m
 mk_zct ='Y';
 sav_fig='N';
 mk_fig ='Y';
 case_1=1;

 if ispc ... % isunix
    hmdir =  [ DRV ':/samhanach/sa01da/'];      
    hmdirM=  [ 'M' ];    
% path_omat=[ 'W:/sa01da/work/minch2/Archive/'];
 else
  hmdir=  ['/home/sa01da/'];      
  hmdirM=  [  hmdir 'M' ];  
 end
 
  addpath([hmdirM '/Mar_Phys/matlab/general/' ]);
  addpath([hmdirM '/Mar_Phys/matlab/seawater/']);
  addpath([hmdirM '/Mar_Phys/matlab/m_map/'   ]); 
  addpath([hmdirM '/Mar_Phys/matlab/t_tide/'  ]); 
  path_omat=[hmdir 'work/minch2/Archive/' ];
  addpath([ '../'  ]); 

   nzl=11;
 % nzl=31;
  data_dir = [this_dir '/test/'];
  
  data_inp=data_dir;
 path_fig=data_inp;
 path_out=data_inp;
 fig_dir =path_fig;

  [status,message,messid] = mkdir (path_out ) ;
  [status,message,messid] = mkdir (path_fig ) ;

   clear status message;
% tic
   caseid1  = 'minch';  modelid1 = '2'; % input model

   caseid  = caseid1 ;modelid = modelid1;

% ## ====================================
        clear cyears iyears;
       mtb=datenum('2016.03.01','yyyy.mm.dd');
       mte=datenum('2017.01.01','yyyy.mm.dd');
   for ii=1:ceil((mte-mtb)/365),...
                    a=datestr(mtb(1)+(ii-1)*365,'yyyy');
       iyears(ii,:)=str2num(a);
   end
   iyears=unique( iyears);
   cyears=num2str(iyears);
% get unique year list;


 case_1=1;         dir_fvm=dir([path_omat '/mat_20*']);
   [nye, ny ]=size(dir_fvm)  ;
          ncy=length(cyears(:,1));

% plt_nest_Dounie ;% with times
                dir_nest=dir([data_dir '/mi*nestnodes.dat']);
 [ks, kc]= size(dir_nest);
ks=1;
for ki=1:ks,...
        clear nest G Gm BC
%    load([ data_dir 'minch2_NL' num2str(ki) '_nest.mat']);
     load([ data_dir 'minch2_Scope'   '_nest.mat']);
     nest.CID = nest.CID_mi;
     nest.mtb = mtb;    nest.mte = mte;
for iy=1:ncy,
  ya=[];  ya=cyears(iy,:);
clear  R OBC dir_mat dir_m
        ny1=1;
 for ny=ny1:nye
      cyear=dir_fvm(ny,1).name(5:8);
      iyear=str2double(cyear);

 if ~strcmp(cyear,ya), continue; end;

     path_M=[path_omat dir_fvm(ny,1).name '/'];

                     dir_mat =dir([path_M '/FV*_all*.mat']);
    [ nf,  lf]= size(dir_mat);
     if nf < 1, continue; end;

%          nf=2
  for lf=1:nf,
         fvnm=[]; fvnm=dir_mat(lf,1).name;
      dir_m(lf,1:length(fvnm))=fvnm;
  end
disp(['ny=' num2str(ny) ' nf=', num2str(nf) ' ; '  dir_fvm(ny,1).name   ]);

for inf=1:nf,...
  clear fl flt;
fl = dir_mat(inf,1).name ;
           flt=fl(11:16); %fl(3:end-4);
     mtime_flt(inf,1)=datenum(flt,'yymmdd');
end
  gof=[];jj=0;
 for ii=1:length(nest.mtb)
  go=[]; go=find(mtime_flt>=nest.mtb(ii) & mtime_flt<=nest.mte(ii));
   if ~isempty(go),
   jj=jj+1;
   gof(jj,1:length(go))=go;
  end
 end
%   if isempty(gof);continue; end
    [ni,nj]=size(gof);
% ## ############################

%% ============= Get FVCOM metrix: 0 ============
% %% ======== ============ obc : 0.1 ========
                       load(['mesh2.mat']);
     Minch2=load([path_M '/' 'mesh.mat']);
     Minch2.mesh.nbe    = mesh.nbe;
     Minch2.mesh.bcells = mesh.bcells;
     Minch2.mesh.bnodes = mesh.bnodes;
 try
     Minch2.mesh.tri   = Minch2.mesh.trinodes;
     Minch2.mesh.x     = Minch2.mesh.nodexy;
     Minch2.mesh.x(:,3)= Minch2.mesh.depth;
 end
 %make z in depth coordinates, not the discrete equivalent of sigman coords.
      nz=nzl-1; %nz=length(FVCOM1.t(1,:,1));
 KBM1=nz;
  Minch2.mesh.nodez=[1:KBM1];
         mesh=Minch2.mesh;
         N= mesh.Nelems ;M=mesh.Nverts;
  dsig=1/(KBM1+1);
 sigvec=(0:(KBM1-1))'*dsig+0.5*dsig;
 for m=1:M
   zall(:,m)=-sigvec*mesh.depth(m);
 end
 for n=1:N
   zuv(:,n)=-sigvec*mesh.uvdepth(n);
 end
 mesh.nodez =zall;
 mesh.zuv   =zuv;
 mesh.sigvec=sigvec;

 Gm.obc_nodes  = nest.NID_mi;
 Gm.obc_elems  = nest.CID_mi;

 Gm.DIM_siglev = length(mesh.nodez(:,1))+1 ; % siglev; 11
 Gm.DIM_siglay = length(mesh.nodez(:,1))+1 ; % siglay; 10
  Gm.zl =[0     1.0 2.0 5.0 10.0      20.0 30.0      50.0      75. 100.       145.];  % 11 % minch
 [Gm.Z,Gm.Z1,Gm.ZZ,Gm.ZZ1,Gm.DZ,Gm.DZ1,Gm.SIGMA] = ...
     setup_sigma_gen0('sigma.dat',Gm.zl,mesh,path_out );
  Gm.D = double(mesh.depth);

   cex=mesh.tri(Gm.obc_elems,:); ncli=length(cex(:,1));
   for ic=1:ncli,
   Gm.xg(ic,1)   = mean(mesh.geog(cex(ic,:),1));
   Gm.xg(ic,2)   = mean(mesh.geog(cex(ic,:),2));
   end

Gm.caseid1  = caseid1 ;
Gm.modelid1 = modelid1;

    G=Gm;
    try nest.CID=nest.CID_mi; end
get_GBC_minch2;
%  [mesh] = cell_area_trans(mesh, nmfcell,I_MFCELL_N);

%======================
%    clear mesh;
% load(['mesh.mat']);  % etive mesh
% get_GBC_etive8;
% G_etive=G;
  BC = [];
%chek if its exists:  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BC_fd=[];
%% ===================
if  case_1, ...
% clear  dir_mat dir_m
      OBC_TIDEOUT_INITIAL =0  ;

[ nf, lf ]= size(dir_m);
  disp(['nf=' num2str(nf) ' mat files to process' ]);
  nmfnode_GL = length(G.obc_nodes);
  nmfcell_GL = length(G.obc_elems);
      KB   = G.DIM_siglev; % 11
      KBM1 = G.DIM_siglay; % KBM1=KB-1;
no=[1:length(G.obc_nodes)]';
ne=[1:length(G.obc_elems)]';

%for jnf=1:nf,... inf=jnf;
 for ii=1:ni,...
  go=[];go=gof(ii,:); go=go(go>0);
  if isempty(go),continue; end

 for jnf=1:length(go),...
  clear R OBC fl flt inf ;
          inf=go(jnf);
  if isempty(inf) | inf <1 , continue; end
                    fl=[];% fl = deblank(dir_m(inf,:));
 					fl = dir_mat(inf,1).name ;
%%  find empty slot:
%            flt=fl(11:16); %fl(3:end-4);
%    mtime_fl=datenum(flt,'yymmdd');
%    mtime_fl=datenum(flt,'yyyymmddHH');
%   datestr(mtime_fl,31)

   disp(['inf=' num2str(inf) ' ... ' fl ]);
%   use chanks

           load([ path_M fl ],'-mat');
  % Roms3d_to_FVcom3d ;% =>R ; OBC
  % Minch2_etive8_3D ; %

  if length(FVCOM1.el(1,:))<=1, FVCOM1.el=FVCOM1.el'; end

    OBC.mtime=FVCOM1.Date;

    a =FVCOM1.el( :,Gm.obc_nodes); OBC.zeta =permute(a,[2,1]) ;
    a =FVCOM1.t(:,:,Gm.obc_nodes); OBC.t    =permute(a,[3,2,1]);
    a =FVCOM1.s(:,:,Gm.obc_nodes); OBC.s    =permute(a,[3,2,1]);
    a =FVCOM1.z(:,:,Gm.obc_nodes); OBC.z    =permute(a,[3,2,1]);
    a =FVCOM1.u(:,:,Gm.obc_elems); OBC.u    =permute(a,[3,2,1]);
    a =FVCOM1.v(:,:,Gm.obc_elems); OBC.v    =permute(a,[3,2,1]);
    a =FVCOM1.ww(:,:,Gm.obc_elems); OBC.ww  =permute(a,[3,2,1]);
    a =FVCOM1.ua(: ,Gm.obc_elems); OBC.rubar=permute(a,[2,1]);
    a =FVCOM1.va(: ,Gm.obc_elems); OBC.rvbar=permute(a,[2,1]);
    a =FVCOM1.h(Gm.obc_nodes);     OBC.h    =a;

 a =FVCOM1.r(          :,:,Gm.obc_nodes); OBC.r             =permute(a,[3,2,1]);
 a =FVCOM1.short_wave(   :,Gm.obc_nodes); OBC.short_wave      =permute(a,[2,1]);
 a =FVCOM1.net_heat_flux(:,Gm.obc_nodes); OBC.net_heat_flux   =permute(a,[2,1]);
 a =FVCOM1.uwind_speed(  :,Gm.obc_elems); OBC.uwind_speed     =permute(a,[2,1]);
 a =FVCOM1.vwind_speed(  :,Gm.obc_elems); OBC.vwind_speed     =permute(a,[2,1]);
 a =FVCOM1.precip(       :,Gm.obc_nodes); OBC.precip          =permute(a,[2,1]);
 a =FVCOM1.evap(         :,Gm.obc_nodes); OBC.evap            =permute(a,[2,1]);

    OBC.ctime=datestr(FVCOM1.Date,31);

       zu=mesh.zuv(:,Gm.obc_elems);
   OBC.zu =permute(zu,[2,1]);

  if  jnf==1,
        BC=OBC;
  else
%        OBC.mtime=BC.mtime+inf;
        jt=length(BC.mtime);
   BC.mtime = [BC.mtime; OBC.mtime];
           ntimes=length(OBC.mtime);

  for it=1:ntimes,...
        jt=jt+1;%it
% BC.time(jt,1)=OBC.mtime(it);
  BC.t(:,:,jt)  =OBC.t(:,:,it);
  BC.s(:,:,jt)  =OBC.s(:,:,it);
  BC.u(:,:,jt)  =OBC.u(:,:,it);
  BC.v(:,:,jt)  =OBC.v(:,:,it);
  BC.ww(:,:,jt) =OBC.ww(:,:,it);
  BC.zeta(:,jt) =OBC.zeta(:,it);
  BC.rubar(:,jt)=OBC.rubar(:,it);
  BC.rvbar(:,jt)=OBC.rvbar(:,it);
  BC.z(  :,:,jt)=OBC.z(:,:,it);

  BC.r(          :,:,jt)= OBC.r(          :,:,it);
  BC.short_wave(   :,jt)= OBC.short_wave(   :,it);
  BC.net_heat_flux(:,jt)= OBC.net_heat_flux(:,it);
  BC.uwind_speed(  :,jt)= OBC.uwind_speed(  :,it);
  BC.vwind_speed(  :,jt)= OBC.vwind_speed(  :,it);
  BC.precip(       :,jt)= OBC.precip(       :,it);
  BC.evap(         :,jt)= OBC.evap(         :,it);

  BC.ctime(jt,:)=OBC.ctime(it,:);
  end
  end
 end; %jnf  in-files


%add one more:
 BC.mtime(end+1)=BC.mtime(end)+mean(diff(BC.mtime(10:15)));
              nt=length(BC.mtime);
              kt=nt-1;
  BC.t(:,:,nt)  =BC.t(:,:,kt);
  BC.s(:,:,nt)  =BC.s(:,:,kt);
  BC.u(:,:,nt)  =BC.u(:,:,kt);
  BC.v(:,:,nt)  =BC.v(:,:,kt);
  BC.ww(:,:,nt) =BC.ww(:,:,kt);
  BC.zeta(:,nt) =BC.zeta(:,kt);
  BC.rubar(:,nt)=BC.rubar(:,kt);
  BC.rvbar(:,nt)=BC.rvbar(:,kt);
  BC.z(  :,:,nt)=BC.z(:,:,kt);

  BC.r(          :,:,nt)= BC.r(          :,:,kt);
  BC.short_wave(   :,nt)= BC.short_wave(   :,kt);
  BC.net_heat_flux(:,nt)= BC.net_heat_flux(:,kt);
  BC.uwind_speed(  :,nt)= BC.uwind_speed(  :,kt);
  BC.vwind_speed(  :,nt)= BC.vwind_speed(  :,kt);
  BC.precip(       :,nt)= BC.precip(       :,kt);
  BC.evap(         :,nt)= BC.evap(         :,kt);

  BC.ctime(nt,:)=datestr(BC.mtime(nt),31);


       OBC_TIDEOUT_INTERVAL=mean( diff(BC.mtime(2:nt-1)) *24*3600) ;
%      OBC_TIDEOUT_INTERVAL=     (BC.mtime(3)-BC.mtime(2))*24*3600 ;
 BC.OBC_TIDEOUT_INTERVAL=OBC_TIDEOUT_INTERVAL;
 BC.OBC_TIDEOUT_INITIAL =OBC_TIDEOUT_INITIAL;

OBC=BC;
%% ===================
  cdates=[datestr(BC.mtime(1),'yyyymmdd') '-' datestr(BC.mtime(end),'yymmddHHMM') ]; %'_' cdates
  save([ path_out 'BC' num2str(ki) '_' cdates '.mat'], 'BC','G','-v7.3' );
% else
%      BC_fd =dir ( [ path_out 'BC_' '*.mat']);
%   load([ path_out BC_fd(end,1).name]);

end; %ii, gof(ii,:)
end; % case 1
end; %ny
end; %iy (from the requested list)
end; % ki
  return

%% ======================1 ===
%   get_zct
