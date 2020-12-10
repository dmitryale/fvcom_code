% function[]=S03_Extract_point_3Duvts_Scope();
% ========================================================================&
% Extracting tidal (1h) FVCOM model currents data  at target location(s)  %
% INput: ['ADCP_Jura..mat';   'minch2_' 'Scope' '_nest.mat'];             &
%        [path_out 'BC1'  '_' cdates '.mat']; % L:6,70                    
% OUTput:[data_dir '/ABC_Kames_JuraWest_2016.mat'], 'FVCOM', 'A'); L145:L318
% [ path_fig '/CMall_FV_20160728-20160815_11.mat' ] extracted u,v,t,s,z,ua
% History : based on R05u_mk_residual_meanflow_1h_NL.m ; 2017.07.21       &
% (c) created : Dmitry Aleynik AT sams.ac.uk, 29/10/2018                  & 
%  __o_O__/X/                                                             &
%  \_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS, COMPASS/CAMPUS 2018-2021     &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

clear all, close all
tic;
% poolSize=1;
%% Loading FVCOM data
this_dir=pwd; DRV=this_dir(1:1);
   iyear=2016; cyear=num2str(iyear);
   
%% Dima's 
 if ispc... % isunix
    hmdir =  [ DRV ':/samhanach/sa01da/'];      
    hmdirM=  [ 'M' ];    
  path_omat=[hmdir '/work/minch2/Archive/'];
% path_omat=[ 'W:/sa01da/work/minch2/Archive/'];

 else
  hmdir=  ['/home/sa01da/'];      
  hmdirM=  [  hmdir 'M' ];  
 end
 
  addpath([hmdirM '/Mar_Phys/matlab/general/' ]);
  addpath([hmdirM '/Mar_Phys/matlab/seawater/']);
% addpath([hmdirM '/Mar_Phys/matlab/m_map/'   ]); 
  addpath([hmdirM '/Mar_Phys/matlab/t_tide/'  ]); 
  path_omat=[hmdir 'work/minch2/Archive/' ];
  addpath([ '../'  ]); 

             data_dir = [this_dir '/'];   %'/Scope/'   
 path_fig = [data_dir '/Figures/'];
[status,message,messageid] = mkdir(path_fig);  clear status message messageid
                path_fvc=[hmdir '/work/minch2/Archive/mat_' cyear '/']; 
 if iyear==2016, path_fvc=[path_omat 'mat_' cyear '_a/'];  end;

%% ========================
global mesh CVM
load(['../mesh.mat'],'-mat');
load([ '../CVM_090101.mat'],'-mat');  
CVM.a1u=double(CVM.a1u); CVM.aw0=double(CVM.aw0); 
CVM.a2u=double(CVM.a2u); CVM.awx=double(CVM.awx); CVM.awy=double(CVM.awy);
try CVM.NV=double(CVM.nv); CVM.NBE=double(CVM.nbe); end
load(['../uk_coastline.mat'],'-mat'); l=find(isnan(coast(:,1)));
%varsg={'t';'s';'ua';'va'}; nv=4; 
% varsg={'t';'s'; 'u';'v';'ww'; 'rubar';'rvbar';'zeta'};   
% nv=length(varsg);
varnames={'mtime';'zeta'; 't'; 's'; 'z'; 'u'; 'v';'ww'; 'rubar'; 'rvbar'; 'r'; ...
'short_wave';'net_heat_flux';'uwind_speed';'vwind_speed'; 'precip'; 'evap';};
  [ nvar,kv]=size(varnames);
% S01_plt_nest_Scope_nodes ; 
      load([ data_dir '/minch2_Scope' '_nest.mat'],'-mat'); % 
% S01b_get_ADCP_data_Scope ;     
      load([ data_dir '/ADCP_Jura_20160728.mat']);  
      
if 1 ...   %L66:L564   
%% Download only that data points from the  FVCOM 
try 
         dir_ma =dir([ data_dir '/BC*' '.mat']);
 load([ data_dir dir_ma(end,1).name ],'-mat');
 load([ data_dir '/minch2_Scope'   '_nest.mat'],'-mat'); % 
try BC.nod = nest.obc_nodes; 
    BC.cel = nest.obc_elems;end
end
%      nest.CID  = nest.CID_mi; 
     ADCP.info =  ADCP.site;
     datestr(ADCP.mtime(1))
  try mesh.sigvec=double(G.Z(1,:)); end % siglays T,S
% try mesh.sigvec=G.ZZ1(1,:); end % siglevs  U,V
     mesh.depth=double(mesh.depth); 
 
%% = get the nodes and lements surrounding single target site
      gxv=ADCP.lon(:)';  % it is fixed !
      gyv=ADCP.lat(:)';
% http://www.nearby.org.uk/coord.cgi?p=114544.6%2C+851118.1&f=full#llNG1454451118;
%     gxv=-5.903137;     gyv=56.064578; 
 
      CM.gxv=gxv';  CM.gyv=gyv';  
    [CMv]=get_xy_elm(CM.gxv(:),CM.gyv(:));
      CM.nod=CMv.nod; 
      CM.elm=CMv.elm; 
    ADCP.cel=CMv.elm;   ADCP.nod = CMv.nod;
 
     ip1=1; xc=mesh.uvnode(:,1); yc=mesh.uvnode(:,2);
%        plot(gxv,gyv,'o'); hold on;
%        plot(xc(ADCP.cel),yc(ADCP.cel),'+'); hold on;
%        plot(mesh.geog(nv(ADCP.cel,:),1),mesh.geog(nv(ADCP.cel,:),2),'v'); hold on;

     [hh] =interp_anodal3Dc(gxv(ip1),gyv(ip1),ADCP.cel(ip1),  mesh.depth, mesh,CVM,xc,yc);  
              CMv.depth=hh;  
     CM.depth=CMv.depth; % zm=mesh.sigvec*hh; layers 
 %   gzv=ADCP.depthi'*CM.depth';     % it is fixed : top, mid-dept, botm       
     
  ibo= [35,26, 1];    zo = ADCP.bin_depth(ibo,1);   %jura
% ibo= [26,17, 3];    zo = ADCP.bin_depth(ibo,1);  %polteit site Skye 
% izm= [ 1, 5,10];    zm = mesh.zuv( izm,ADCP.cel);% zm=[4.28,  12.84, 27.10];
            % it is fixed : top, mid-dept, botm       
             gzv=zo;
  %get all levels
  izm= [ 1:10];% zm = mesh.zuv( izm,ADCP.cel);% zm=[4.28,  12.84, 27.10];
                 zm = mesh.sigvec*hh;
         gzv=abs(zm(:)); % all layers
  CM.gzv    =gzv;
  ADCP.gzv  =gzv;
  ADCP.depth=gzv;%NB
  try  ADCP.depthi= ADCP.bin_depth(:,1); end
 
  nbe=CVM.NBE;
   nv=CVM.NV;
         tri_ids=CMv.elm;
      ii=tri_ids;
% Surrounding Element IDs
    e1  = nbe(ii,1);
    e2  = nbe(ii,2);
    e3  = nbe(ii,3);      
% !element location (i) and surrounding nodes (n1,n2,n3)
     n1  = nv(ii,1);
     n2  = nv(ii,2);
     n3  = nv(ii,3);
CM.nodes=[n1,n2,n3];
CM.elems=[e1,e2,e3];
nodes=[CM.nodes,CM.nod];
elems=[CM.elems,CM.elm]; elems=elems(elems>0); CM.elems=elems;
%%  ================

%  ADCP.mtime=BC.mtime;
%% need sigloc and metrics; zz1 dz1
try mesh.tri   =mesh.trinodes; end
try mesh.x     =mesh.geog;  end
try mesh.x(:,3)=mesh.depth; end

get_siglocs; % setup_sigma_gen;
if 0....         
 if 1, %L145 get the ABC data from already extracted nearest points
  %% Download only that data points from the  FVCOM 
     FVCOM.Date = BC.mtime;
%% select only those who fit the CM metrix !!!
    elemsu=elems(:);    nodesu=nodes(:);
    clear  nodBC elmBC 
    k=0;
for ii=1:length(nodesu),...(elemsu)
    go=[]; go= find (nodesu(ii) == BC.nod);    
  if ~isempty(go), k=k+1;
  nodBC(k,1) =go;
  end
end
    k=0;
for ii=1:length(elemsu),...    
    go=[]; go = find (elemsu(ii) == BC.cel);
  if ~isempty(go),k=k+1;
    elmBC(k,1) = go ;
  end
end

 FVCOM0.nodes =BC.nod(nodBC);
 FVCOM0.elems =BC.cel(elmBC);
 a= BC.zeta(nodBC,:); FVCOM0.el=a'; 
 a = BC.u(elmBC,:,:); FVCOM0.u =permute(a,[3,2,1]);
 a = BC.v(elmBC,:,:); FVCOM0.v =permute(a,[3,2,1]);
 try a = BC.ww(elmBC,:,:); FVCOM0.ww =permute(a,[3,2,1]); end
 a = BC.rubar(elmBC,:,:); FVCOM0.ua =permute(a,[2,1]);
 a = BC.rvbar(elmBC,:,:); FVCOM0.va =permute(a,[2,1]);
 a = BC.t(nodBC,:,:); FVCOM0.t =permute(a,[3,2,1]);
 a = BC.s(nodBC,:,:); FVCOM0.s =permute(a,[3,2,1]);
 a = BC.uwind_speed(elmBC,:); FVCOM0.uwind_speed =permute(a,[2,1]);
 a = BC.vwind_speed(elmBC,:); FVCOM0.vwind_speed =permute(a,[2,1]);
     a= BC.h(nodBC); FVCOM0.hh=a'; 
 FVCOM0.Date=BC.mtime;
     FVCOM=FVCOM0;     
 else % L145:L181
%%% ========================== II all fvcom files ============   
   try ADCP.mtime=BC.mtime; end
       fv_mt=ADCP.mtime; 
ibf=find(fv_mt >= min(ADCP.mtime)-1.0,1,'first');
ief=find(fv_mt <= max(ADCP.mtime)+1.25,1,'last' );
 fv_go  =[ibf:ief]; fv_go1 = fv_go;
%% extract the first file only for getting metrics 
   cy=datestr(fv_mt(1),'yyyy');
  
    dir_fvm=dir( [ path_omat '/mat_20*']);
%   dir_fvc=dir( [ path_omat dir_fvm(end).name '/FV*.mat' ]);
   for jy=1:length(dir_fvm),...
        cyf=dir_fvm(jy).name(5:8);   
    if strcmpi(cy, cyf), iy=jy; end
   end
   
    dir_fvc=dir( [ path_omat dir_fvm(iy).name '/FV*.mat' ]);
    
           fnm=[ dir_fvm(iy).name '/' dir_fvc(end).name ];
    A{1,1}=load([ path_omat fnm]);   
%% =======================
FVCOM=A{1,1}.FVCOM1; Date=FVCOM.Date; 
  
   av=fieldnames(A{1,1}.FVCOM1); bv=av(19:length(av))'; 
                bv=cat(2,bv,{'time';'Itime';'Itime2';}'); % 'r';'z'
FVCOM=rmfield(FVCOM,bv);
vars_fv={'t';'s';'u';'v';'ww';'ua';'va';'el';}; nva=length(vars_fv);%'r';'z';;;'ww'
   na=length(A);%  length(dir_fvc); 
if na>1,...
 for ki=2:na
    b=A{ki,1}.FVCOM1.Date;
    go=find(b>FVCOM.Date(end));
    if isempty(go),continue;end
    FVCOM.Date=cat(1,FVCOM.Date,b(go));
    for vi=1:nva
        va=char(vars_fv(vi));
  if ~isfield (FVCOM,[va]), disp([ ' ' va ' : do not exists in FVCOM L188']); continue; end
        
        a= FVCOM.(va);
        b= A{ki,1}.FVCOM1.(va);
        c= size(a);
        if length(c)==3, a=cat(1,a,b); end
        if length(c)==2, a=cat(1,a,b); end
        FVCOM.(va)=a;
    end
 end
end

path_out=path_fig;%  = ['../Figures/'];
Date=FVCOM.Date;
    [ri,rc]=size(ADCP.lon);
     go_A=[1:length(ADCP.mtime)]'; n=length(go_A);      
        if n>rc,
            go_A=[2:n-1]';
        end
        
mt_obs=[ floor(fv_mt(1)) (ceil(fv_mt(end)))] ;datestr (mt_obs)
 
%%=====================
%% extract all files now for getting ALL FVCOM data
%  removing unnesessary bits
   
for ki=1:length(dir_fvc(:))...% fv_go1),...
                ik=fv_go1(ki); 
%    fnm=dir_fvc(ik).name; dir_fvm(end).name'/'
     fnm=[  dir_fvc(ki).name ];
     mt_fv(ki,1) = datenum(fnm(11:16),'yymmdd');
end
  good_fv=[]; good_fv=find(mt_fv >= mt_obs(1) & mt_fv <= mt_obs(end) );
      
    kj=0;
for ki=1:  length(good_fv)  ,...  
           kj=good_fv(ki);
    fnm=[  dir_fvc(kj).name ];
    A{ki,1}=load([path_fvc fnm]);    
    A{ki,1}.FVCOM1=rmfield( A{ki,1}.FVCOM1,bv);
    A{ki,1}.FVCOM1.filename=fnm;
end

%% =======================
FVCOM=A{1,1}.FVCOM1; Date=FVCOM.Date; 
na=length(A); % av=fieldnames(A{1,1}.FVCOM1); bv=av(19:length(av)); 'ww'
 try FVCOM=rmfield(FVCOM,{'time';'Itime';'Itime2';'r';'z';}); end; 
 try FVCOM=rmfield(FVCOM,bv); end
   FVCOM1=FVCOM;
            FVCOM=rmfield(FVCOM,vars_fv);            
            
%     nodes=ADCP.nod;        elems=ADCP.cel;         %single point
    nodes=nest.obc_nodes;  elems=nest.obc_elems;  % six elms, 5 nodes

 FVCOM.nodes = nodes;          
 FVCOM.elems = elems;        
    fnm=A{1,1}.FVCOM1.filename;
   try r1=find(fnm=='/',1,'last')+1; end
    r1=find(fnm=='_',1,'first')+4;
    r2=find(fnm=='.',1,'last')-1;
      fnms=A{1,1}.FVCOM1.filename(r1:r2);
FVCOM.fnm(1,1:length(fnms))=fnms;

    for vi=1:nva,
        va=char(vars_fv(vi));    pnts=[];    
   if ~isfield (FVCOM1,[va]), disp([ ' ' va ' : do not exists in FVCOM1 L218']); continue; end
        if (vi> 2 | vi < 8 ), pnts=elems; end
        if (vi<=2 | vi > 7 ), pnts=nodes; end
%         if vi==nva, stop; end
          a = FVCOM1.(va);
        c= size(a);
        if length(c)==3, b=squeeze(a(:,:,pnts)); end
        if length(c)==2, b=squeeze(a(:,pnts))  ; end
        FVCOM.(va)=b;    
    end

if na>1,...
        
 for ki=1:na,...
            b=A{ki,1}.FVCOM1.Date;
    go=find(b>FVCOM.Date(end));
    if isempty(go),continue;end
                fnm=A{ki,1}.FVCOM1.filename(r1:r2);
FVCOM.fnm(ki,:)=fnm;    
    FVCOM.Date=cat(1,FVCOM.Date,b(go));    
    for vi=1:nva
         va=char(vars_fv(vi));    pnts=[];    
   if ~isfield (FVCOM1,[va]), disp([ ' ' va ' : do not exists in FVCOM1 L218']); continue; end
        if (vi> 2 | vi < 8 ), pnts=elems; end
        if (vi<=2 | vi > 7 ), pnts=nodes; end
          a= FVCOM.(va);
        b= A{ki,1}.FVCOM1.(va);
        c= size(b);
        if length(c)==3, a=cat(1,a,squeeze(b(:,:,pnts))); end
        if length(c)==2, a=cat(1,a,squeeze(b(:  ,pnts))); end
        FVCOM.(va)=a;
    end
 end
end
     save([data_dir '/ABC_Kames_JuraWest_2016.mat'], 'FVCOM', 'A'); % ,'A','CM'  
 end; %if L145:318
else
    
   disp([' continue working with old ABC file and new CM.elm: ' num2str(CM.elm) ,' '])
end
 %% ====================
     get_siglocs; % setup_sigma_gen; % updATE IT
 ri=length(ADCP.depth); 
 kis=1:ri;  % number of rows (levels in ADCP);
  IZ=1;     % consider 1 level only at first instance
  vars_fv={'t';'s';'u';'v';'ww';'ua';'va';'el';}; nva=length(vars_fv);%'r';'z';;;'ww'

   clear go mt mtv a b c;
%  clear A ;
     xc=mesh.uvnode (:,1);
     yc=mesh.uvnode (:,2);

                  Date=FVCOM.Date;
   [nt,nz,N]=size(FVCOM.u); lvls=nz ; % kbm1
%%
      go_A=[1:length(ADCP.mtime)]'; n=length(go_A);      
    ip1=1;
% plot_grid_pnt 
 FVCOM.xxn(1,:)=mesh.geog(FVCOM.nodes(:),1);   FVCOM.yyn(1,:)=mesh.geog(FVCOM.nodes(:),2);
 FVCOM.xxc(1,:)=mesh.uvnode(FVCOM.elems(:),1); FVCOM.yyc(1,:)=mesh.uvnode(FVCOM.elems(:),2);  
%  plot(FVCOM.xxn,FVCOM.yyn,'-s') ;hold on;
%  plot(FVCOM.xxc,FVCOM.yyc,'^') ;hold on;
%  plot(ADCP.lon,ADCP.lat,'+') ;hold on;
 
     clear CMal;
for ip1=1:length(gxv),...     
     CMal(ip1).CM=CM;        
     clear hl hlg gl_id  dpth cdpth 
   
for ki=kis, 
    IZ=ki; 
    clear gl_id dpth cdpth       
    gl_id=[deblank(ADCP.info(ip1,1:8)) '_iz_' num2str(IZ,'%2.2d') ];
% UVa-------------
    mtv=ADCP.mtime(go_A);            % nt=length(mtv);
          dto=mean(diff(mtv(2:10)));  % dto~=1/24;   % 0.3333/24 ? obs
          dtm=mean(diff(Date(2:10))); % dtm=1/24;    
    go_g=find(mtv >=Date(1) & mtv<=Date(end)); %obs
    go_f=find(Date>=mtv(1)  & Date<=mtv(end)); %model
    disp(['UVa: '  num2str(ki) ' go_g=' num2str(length(go_g)) ' go_f=' num2str(length(go_f))]);
    if isempty(go_g) | isempty(go_f), continue; end
    
    CM.go_fv=go_f;
    CM.go_gv=go_g;
    CM.fi_mtv=Date(go_f);
    CM.g_mtv = mtv(go_g);

%   vib=3; vie=4; 
    vib=3;vie=vib+0;  % u only & zeta
    for vi=vib:vie,...
        va=char(vars_fv(vi));
  if ~isfield (FVCOM,[va]), disp([ ' ' va ' : do not exists in FVCOM  L284']); continue; end
          a=FVCOM.(va);           % 2D or 3D
%         b=ADCP.([va])(IZ,go_g); % 1D
        tri_ids=CM.elm;
   [hh] =interp_anodal3Dc(gxv(ip1),gyv(ip1),tri_ids(ip1),       mesh.depth, mesh,CVM,xc,yc); %2D
    CM.(['fi_h_' va(1)])(ip1,IZ)=hh ;   
                 npv=length(mtv);   
            if ispc, 
               npv=length(go_g); %npv=49; 
            end;                 
        for ip=1:npv,...                
                %  minimize time for that ip==jt
            clear jt mtvi go_fi vel vel_1 vel_2 vel_3  el1 el2 el3
            jt=go_g(ip);
          mtvi= mtv(jt);
          
               go_fi=find( Date >= mtvi-dtm & Date < mtvi+dtm  );
%              go_fi=find( Date >= mtvi-dto & Date < mtvi+dto  );
            if isempty (go_fi) || length(go_fi)<1, continue; end
           
    if vi==vib ,   
            factm_a=0; factm_b=0;
            factm_a = 24*(mtvi-Date(go_fi(1))); 
            if length(go_fi)==1             
               elev=zeros(1,mesh.Nverts); 
               elev(:,FVCOM.nodes)=FVCOM.el(go_fi,:);                
                [el3]=interp_anodal3Dc(gxv(ip1),gyv(ip1),tri_ids(ip1),squeeze(elev(1,:)),mesh,CVM,xc,yc); %2D
              else                  
               factm_b = 24*(Date(go_fi(2))-mtvi);
               elev=zeros(2,mesh.Nverts); 
               elev(:,FVCOM.nodes)=FVCOM.el(go_fi,:);   
% test the location :               
% plot (     gxv(ip1),gyv(ip1),'d'); hold on;
% set(gca,'xlim',[gxv(ip1)-0.01 gxv(ip1)+0.01] );  set(gca,'ylim',[gyv(ip1)-0.01 gyv(ip1)+0.01]);
% plot ( mesh.uvnode(tri_ids(ip1),1), mesh.uvnode(tri_ids(ip1),2),'r^'); hold on;
% plot ( mesh.geog(mesh.trinodes(tri_ids(ip1),:),1), mesh.geog(mesh.trinodes(tri_ids(ip1),:),2),'gs'); hold on;
% jnk=patch('Vertices',[mesh.geog(:,1),mesh.geog(:,2)],'Faces',mesh.tri,...
% 	      'EdgeColor',0.7*[1 1 1];,'FaceColor','none'); hold on;  

               [el1]=interp_anodal3Dc(gxv(ip1),gyv(ip1),tri_ids(ip1),squeeze(elev(1,:)),mesh,CVM,xc,yc); %2D
               [el2]=interp_anodal3Dc(gxv(ip1),gyv(ip1),tri_ids(ip1),squeeze(elev(2,:)),mesh,CVM,xc,yc); %2D
                el3=el1*(1-factm_a) + el2*(1-factm_b);
               end  % um=zeros(2,mesh.Nelems);             
                CM.(['fi_el' va(1)])(ip,IZ)=el3 ;
                nzr=1;
                z_tmp= -mesh.siglev*(hh+el3);
%               sig_loc=interp1([0;z_tmp],[0;mesh.siglev],gzv(ip));
                   sig_loc=interp1([z_tmp],[mesh.siglev],gzv(IZ,ip1));
                if sig_loc==0,      sig_loc=mesh.siglev(1);   end
                if isnan(sig_loc),  sig_loc=mesh.siglev(end) ;end
                CM.(['fi_sig_loc' va(1)])(ip,IZ)=sig_loc;
            else
                sig_loc = CM.('fi_sig_locu')(ip,IZ);
            end
            
        end; %np
    end; %vi
    
%%  apply to U,V but skip T,S:================
%% {{{{{{{{{{{{{{{{{{{{{{{{{{
%     mt=ADCP.mtime;  
      mt=ADCP.mtime(go_A);  ntA=length(mt);
    go_g=find(mt>=Date(1) & mt <= Date(end)); % 20 minutes
    go_f=find(Date>=mt(1) & Date <= mt(end)); 
    CM.gl_id=gl_id;           
    disp( ['U,V : ' num2str(ki) ' go_g=' num2str(length(go_g)) ' go_f=' num2str(length(go_f)) ' ' gl_id ]);
        if isempty(go_g) | isempty(go_f), continue; end
    CM.go_f=go_f; CM.go_g=go_g;
    CM.fi_mt=Date(go_f);
    CM.g_mt=mt(go_g);
    dba=datestr(CM.g_mtv( 1),'yyyymmdd_HHMM');
%     save([ path_out '/' gl_id '_FV_u' dba  '.mat'],'CM','-v7.3');
    
    %% ===========================================
    gx=gxv(ip1); % ADCP.lon(ip1)';
    gy=gyv(ip1); % ADCP.lat(ip1)';
    gz=ADCP.depth(IZ,ip1)';
    
    CM.gx(ip1)=gx; CM.gy(ip1)=gy; CM.gz(ip1,IZ)=gz;
    dtm=mean(diff(Date(2:10))); % dtm=1/24;    
    
                   npv=length(go_g);
%       if ispc,   npv=49; end;
    %% =================     vib=1;vie=2; %t,s
    vib=3;vie=vib+1;  % u,v
  for vi=vib:vie,...
      va=char(vars_fv(vi));
  if ~isfield (FVCOM,[va]), disp([ ' ' va ' : do not exists in FVCOM L353']); continue; end
      a = FVCOM.(va);          
%     b = ADCP.(va)(IZ,go_g)';         
%     CM.(['g_' va])(:,IZ) = b;
    for ip=1: npv,...           
            clear jt mtvi go_fi vel vel_1 vel_2 vel_3 factm_a factm_b 
            jt=go_g(ip);
          mtvi= mt(jt);          
                 go_fi=find( Date >= mtvi-dtm & Date < mtvi+dtm   );                
    if isempty (go_fi) || length(go_fi)<1, continue; end
    
            factm_a=0; factm_b=0;
            factm_a = 24*(mtvi-Date(go_fi(1))); 
            sig_loc = CM.('fi_sig_locu')(ip,IZ);
         if sig_loc==0, sig_loc=mesh.siglev(1)  ; end;
   if isnan(sig_loc),   sig_loc=mesh.siglev(end); end;
            if length(go_fi)==1             
      val=zeros(1,nz,mesh.Nelems); 
      val(:,:,FVCOM.elems)=double(squeeze(a(go_fi,:,:))); %3D u,v     
             aa=squeeze(val(1,:,:))';      
%    [val_3]=interp_azonal3D(gxv(ip1),gyv(ip1),sig_loc,lvls,tri_ids(ip1),aa,mesh, CVM  ); %u,v           
     [val_3]=interp_azonal3Dc(gxv(ip1),gyv(ip1),sig_loc,lvls,tri_ids(ip1),aa,mesh, CVM,xc,yc  ); %u,v           
            else                                    
               factm_b = 24*(Date(go_fi(2))-mtvi);
      val=zeros(2,nz,mesh.Nelems);            
      val(:,:,FVCOM.elems)=double(squeeze(a(go_fi,:,:))); %3D u,v
      aa=squeeze(val(1,:,:))';
      bb=squeeze(val(2,:,:))';
%              varargin={gxv(ip1),gyv(ip1),sig_loc,lvls,tri_ids(ip1),aa,mesh, CVM};
 [val_1]=interp_azonal3Dc(gxv(ip1),gyv(ip1),sig_loc,lvls,tri_ids(ip1),aa,mesh, CVM,xc,yc  ); %u,v
 [val_2]=interp_azonal3Dc(gxv(ip1),gyv(ip1),sig_loc,lvls,tri_ids(ip1),bb,mesh, CVM,xc,yc  );  
         val_3 = val_1*(1-factm_a) + val_2*(1-factm_b);
            end
    CM.(['fi_' va])(ip,IZ)=val_3 ;
   end;%ip
  end; %vi

%% }}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
      npv=length(CM.fi_u(:,IZ));
      ciz = num2str(IZ,'%2.2d');
     dpth = round(mean(ADCP.depth(IZ,ip1)));
    cdpth = num2str(dpth,'%3.3d');          
    dba=datestr(CM.g_mt(1),'yyyymmdd_HHMM');
    dbe=datestr(CM.g_mt(npv),'yyyymmdd_HHMM');
                                    dbae=[dba(1:8) '-' dbe(1:8)];
    save([path_fig '/' gl_id '_FV_' dbae  '.mat'],'CM','-v7.3');
  if  1,    
    xlmt=[min(CM.g_mt) max(CM.g_mt) ];
    if diff(xlmt)>3, xlmt=[xlmt(1) xlmt(1)+2.5]; end
    a=[];% a=cat(1,a,CM.g_u(:,IZ) ); a=cat(1,a,CM.g_v(:,IZ)); 
          a=cat(1,a,CM.fi_u(:,IZ)); a=cat(1,a,CM.fi_v(:,IZ));
          aa=max(abs(a)); 
              ylmta=double([-aa aa ]);
    if IZ==1, ylmt=ylmta;end
        ylmt=[-max(max(abs(ylmt),aa)),max(max(abs(ylmt),aa)) ];
        ylmt=double(ylmt);
    if ki==1, figure(ip1+100);clf; end
    
     colr=['brgcmyk'];   colr=cat(1,colr,colr);
     col=colr(ki);
    subplot(2,1,1);
%   plot(CM.g_mt, CM.g_u(:,IZ),'ro-','markersize',2);hold on;
    plot(CM.g_mt(1:npv), CM.fi_u(:,IZ),[col 's-'],'markersize',1.5);hold on;
    set(gca,'ylim',ylmt); set(gca,'xlim',xlmt);    grid on;
    glid=[]; glid=gl_id; glid(glid=='_')='-';
    cva=char(vars_fv(3)); ylabel(['Obs. ' glid ', ' cva ', m/s' ],'color','r');
    datetick('x','keeplimits');
    set(gca,'YMinorTick','on','XMinorTick','on');
 
    subplot(2,1,2);
%   plot(CM.g_mt, CM.g_v(:,IZ),'ro-','markersize',2);hold on;
 hc(ki)=plot(CM.g_mt(1:npv),CM.fi_v(:,IZ), [col 's-'],'markersize',1.5);hold on;
     set(gca,'ylim',ylmt); set(gca,'xlim',xlmt);   grid on; 
                         cva=char(vars_fv(4));
     ylabel([ 'FVcom: ' cdpth ' m, ' cva ' , m/s' ],'color','b');
    datetick('x','keeplimits');
  set(gca,'YMinorTick','on','XMinorTick','on');
                 a=[];  a=[' ' num2str(IZ,'%2.2d') ', Z=' cdpth,'m' ] ;
  hl(ki,1:length(a))=a;
    if ki==kis(end);
        hlg = legend (hc, {hl});
        set(hlg,'location','best'); set(hlg,'box','off');
        set(hlg,'position',[ 0.65 0.40 0.24 0.13]);        
        text(xlmt(1)+0.1,ylmt(1)+0.15,dba(1:8));
    end
 set(gcf,'color','w'); % 'none')
 orient portrait ;
 set(gcf,'PaperPositionMode','auto');
 set(gcf,'renderer','zbuffer') ;
     fig_names  =([ path_fig '/' gl_id '_FV_' dbae  ]);
  print(['-f'],'-dpng','-loose','-r300',[fig_names '.png']);
  end
%   figure(22);clf 
%   a=CM.fi_u(:,1);b=CM.fi_v(:,1); scatter(a,b); axis equal
 
end;%ki
     CMal(ip1).CM=CM;        
end; %ip1
 
    save([path_fig '/CMall_FV_' dbae '_' num2str(length(CM.gzv)) '.mat'],'CMal','-v7.3');
toc;
%   pcolor(FVCOM.u(:,:,4)'); shading interp;  colormap(vivid('br')); colorbar
%   pcolor(FVCOM.u(:,:,4)'); shading interp;  colormap(flipud(vivid('br'))); colorbar
end