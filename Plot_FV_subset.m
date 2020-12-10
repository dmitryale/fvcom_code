%function[]=plot_FV_subset();
%(c) dmitry.aleynik@sams.ac.uk 2017.09.01                                 & 
%  __o_O__¬                                                    FASTNEt    &
%  /_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS Glider mission 4, Malin Shelf &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 clear all 
 close all
  this_dir=pwd; DRV=this_dir(1:1);
     Tb=datenum('2017-09-06 08:00:00',31); % dye release 1kg rodamin 
     Te=datenum('2017-09-07 16:00:00',31);
tic;
addpath M:/Mar_Phys/matlab/m_map
addpath M:/Mar_Phys/matlab/seawater
addpath M:/Mar_Phys/matlab/general/
        a=load('etive6_0coast.mat');
 path_fig='../fig/';
 path_mat='../mat/';

 filePO=['Poltips_Oban_20170901-0930.txt'];
 [PO] = import_pol_oban(filePO);
  v_po = datevec(PO.mtime(1)   );
    time_offset = 678942; % from FVCOM time to matlab time
   [IM_range_lat ,IM_range_lon] = deal([56.445 56.496667 ],[ -(5+28.0/60) -(5+22.5/60)]); cm='_AM'; %Ardmucknish
   ml2m=1852*60; rlat= mean(IM_range_lat);
   xv=[IM_range_lon];   yv=[IM_range_lat];
   m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','wgs84') 
 load('A_Grid.mat','-mat');
 %% ======================
 %% Get gridded subsets                
if 1,...
 now1=now;   disp([ ' Loading Start at: '  datestr(now1,31)]);
  dir_fvmg=dir ([path_mat 'FV_Grd20*.mat']);
%  vars={'el';'t';'s'; 'u';'v';'ww';'ua';'va';'uwind_speed';'vwind_speed'};
   vars={'el';'u';'v';'ww';}; 
 [fn, fi]=size(dir_fvmg);
%   fn=10; % 30 h, the latest files 11-13 are slow to get at laptop
%   fn=13; % 32 h, 48.230603 seconds  at sams PC
for fi=1:fn,...
   A=[];got=[];
   A=load([path_mat dir_fvmg(fi).name ]);
   got=[];  got= find(A.FV.mt >=Tb & A.FV.mt <=Te);
   if isempty(got), continue;end
 %% define vertical levels etc :  
 if fi==1,...
  FV=[];  mt=[];
  FV.mask=A.FV.mask; FV.x=A.FV.x; FV.y=A.FV.y; FV.siglev=A.FV.siglev;
  FV.h=A.FV.h;         mt=A.FV.mt;
% [      nx,ny] = size(A.FV.h);
% [nt,   nx,ny] = size(A.FV.el);
  [nt,nz,nx,ny] = size(A.FV.u);
          nt = length(got);
       FV.mt=A.FV.mt;
       FV.el=A.FV.el;
       FV.u=A.FV.u;
       FV.v=A.FV.v;
       FV.ww=A.FV.ww;
       FV.fi(fi,1)=fi;
       FV.fnm(fi,:)=dir_fvmg(fi).name ;
       FV.go_nt(fi,1)=length(got);
 else
    go=[];  go = find(A.FV.mt > FV.mt(end));
       FV.mt=cat(1,FV.mt,A.FV.mt(go)     );
       FV.el=cat(1,FV.el,A.FV.el(go,:,:)  );
       FV.u =cat(1,FV.u ,A.FV.u( go,:,:,:)); 
       FV.v =cat(1,FV.v ,A.FV.v( go,:,:,:)); 
       FV.ww=cat(1,FV.ww,A.FV.ww(go,:,:,:)); 
       FV.fi(fi,1)=fi;
       FV.go_nt(fi,1)=length(go);
       FV.fnm(fi,:)=dir_fvmg(fi).name ;
 end 
 disp([ ' Loading ' num2str(fi) ' ' dir_fvmg(fi).name ' at: '  datestr(now,31),' dt=' num2str((now-now1)*24*3600) ' s']); 
end; %fn
toc ;   % ~100 seconds over 10 files , 30 muntes for 13 files 
disp([ ' Loading End at: '  datestr(now,31),' dt=' num2str((now-now1)*24*3600) ' sec']); 
clear A; 

%% get the vertical sigma-layers depth : zs  
    nt=length(FV.mt);   % tic
   zs=zeros(nt,nz+1,nxny);
for it=1:nt,...,
 for i=1:nx,
 for j=1:ny,...
   if FV.mask(i,j)==0, continue; end
       yz = -( -( FV.h(i,j)+FV.el(it,i,j) ).*FV.siglev(1:end) ) + FV.el(it,i,j);     
   zs(it,:,i,j) =yz;
 end; 
 end;
end
   FV.zs=zs; 
   clear zs; % toc % 56 sec;
end

 plt_FV_subset
 
%  plt_FV_subset_vect
