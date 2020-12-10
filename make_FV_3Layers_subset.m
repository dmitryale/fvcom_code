%function[]=make_FV_3Layers_subset();
%(c) dmitry.aleynik@sams.ac.uk 2017.09.01                                 & 
%  __o_O__¬                                                    FASTNEt    &
%  /_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS Glider mission 4, Malin Shelf &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
if 1,...
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
    [nz ~]=size(FV.siglev); nz=nz-1;
  yz=[]; zs=[]; zs=zeros(nt,nz+1,nx,ny);
for it=1:nt,...
 for i=1:nx,
 for j=1:ny,...
   if FV.mask(i,j)==0, continue; end
%  yz=zeros(nz+1,1);
   yz=-( -( FV.h(i,j)+FV.el(it,i,j) ).*FV.siglev(1:end ) ) + FV.el(it,i,j);
 %     -( -(  36       +2            ) * [0:-1]           ) + 2 
   zs(it,:,i,j) =yz;
 end; 
 end;
end
   FV.zs=zs; 
   clear zs; % toc % 56 sec;
 disp([ ' Zlayers End at: '  datestr(now,31),' dt=' num2str((now-now1)*24*3600) ' sec']); 
   
end;

FV3=FV;
FV3=rmfield(FV3,{'u';'v';'ww';'zs'});
FV3.layers=[1,fix(nz/2),nz];
for iz=1:length(FV3.layers),...    
        jz=FV3.layers(iz);
%   a=[];  a=FV.zs(:,jz,:,:);    
FV3.zl(:,iz,:,:)=FV.zs(:,jz,:,:);  ... end;      
FV3.u( :,iz,:,:)=FV.u(:,jz,:,:);
FV3.v( :,iz,:,:)=FV.v(:,jz,:,:);
FV3.ww( :,iz,:,:)=FV.ww(:,jz,:,:);
end


  dir_fvmg=dir ([path_mat 'FV_Grd20*.mat']);
  db=datestr(FV3.mt(1),'yyyymmdd-HHMM');
  de=datestr(FV3.mt(end),'yyyymmdd-HHMM');
save([path_mat 'FV_3Layers_Grd_' db '-' de  '.mat'],'FV3','-v7.3');

disp([ ' 3 Layers saved at: '  datestr(now,31),' dt=' num2str((now-now1)*24*3600) ' sec']); 
% test it here 

  a=FV3.zl(:,:,10,10);
   b=FV.zs(:,:,10,10);
  plot(FV3.mt,b); hold on;
  plot(FV3.mt,a,'linewidth',2); hold on;
  plot(FV3.mt,FV3.el(:,10,10),'k'); hold on;
  set(gca,'xlim',[min(FV.mt) max(FV.mt)]);
  ylabel('depth,m');
   set(gca,'YMinorTick','on','XMinorTick','on');

  datetick('x','keeplimits');
   print(['-f'],'-dpng','-loose','-r300',['z_sigma_i10j10'  '.png']);
  
 if 0 
 %% remake FV-3hours files  got=[];  got= find(FV.mt >=Tb & FV.mt <=Te);
      nt=length(got);mt=FV.mt;
 got0 =  find(mod(mt(:),3/24)==0);
 got0(end+1)=got(end);
        nk=length(got0);
 FV0=FV; vars={'el';'u';'v';'ww';};

 for ik=1:nk,...
     if ik==1,
         ib=1;
     else
         ib=got0(ik-1)+1;
     end
         ie=got0(ik);
       go=[ib:ie]'; ng =length(go);
       FV=[];
       FV=FV0;
       FV=rmfield(FV,{'u';'v';'ww';'zs';'el',;'mt'});
FV.mt(1:ng,1)    =FV0.mt(go,1)   ;    
FV.el(1:ng,:,:  )=FV0.el(go,:,:  );    
FV.u( 1:ng,:,:,:)=FV0.u( go,:,:,:);
FV.v( 1:ng,:,:,:)=FV0.v( go,:,:,:);
FV.ww(1:ng,:,:,:)=FV0.ww(go,:,:,:);
% FV.zs(1:ng,:,:,:)=FV0.zs(go,:,:,:);       

    db=datestr(FV.mt(1),'yyyymmdd-HHMM');
    de=datestr(FV.mt(end),'yyyymmdd-HHMM');
    FVgrd_nm=['FV_Grd' db '_' de];
    save([path_mat   FVgrd_nm '.mat'],'FV','-v7.3');
    disp(['saved :'  FVgrd_nm ', nt=', num2str(length(FV.mt)),'; ' , datestr(now,31)]);
 end
end
 
  
  %  plt_FV_subset
 
%  plt_FV_subset_vect
