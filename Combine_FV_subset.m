%function[]=Combine_FV_subset();
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
 path_matg='../mat/';
 path_mat=([DRV ':\HECTOR\matlab\etive8_results\plots\etive27_2013_e\']);

 filePO=['Poltips_Oban_20170901-0930.txt'];
 [PO] = import_pol_oban(filePO);
  v_po = datevec(PO.mtime(1)   );
    time_offset = 678942; % from FVCOM time to matlab time
% [IM_range_lat ,IM_range_lon] = deal([56.433 56.50 ],[ -(5+28.0/60) -(5+22.5/60)]); cm='_AM'; %Ardmucknish
  [IM_range_lat ,IM_range_lon] = deal([56.445 56.496667 ],[ -(5+28.0/60) -(5+22.5/60)]); cm='_AM'; %Ardmucknish
   ml2m=1852*60; rlat= mean(IM_range_lat);
   xv=[IM_range_lon];   yv=[IM_range_lat];
   m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','wgs84') 
 load('A_Grid.mat','-mat');
 %% ======================
 
now1=now;
  disp([ 'start at: '  datestr(now1,31)]);
%%gridded subsets
                 path_mato=path_fig;
  dir_fvmg=dir ([path_mato 'FV_Grd20*.mat']);
  
%  vars={'el';'t';'s'; 'u';'v';'ww';'ua';'va';'uwind_speed';'vwind_speed'};
   vars={'el';'u';'v';'ww';};
   
 [fn, fi]=size(dir_fvmg);
    k=0;
for fi=1:fn,...
   A=load([path_mato dir_fvmg(fi).name ]);
   got=[];  got= find(A.FV.mt >=Tb & A.FV.mt <=Te);
   if isempty(got), continue;end
 
%% define vertical levels:
   FV=[];  mt=[];
% FV.mask=G.mask;   FV.x=G.xq;   FV.y=G.yq;   FV.siglev=G.siglev;
 FV.mask=A.FV.mask; FV.x=A.FV.x; FV.y=A.FV.y; FV.siglev=A.FV.siglev;
 FV.h=A.FV.h;          mt=A.FV.mt;
% [      nx,ny ]=size(A.FV.h);
% [nt,   nx,ny ]=size(A.FV.el);
  [nt,nz,nx,ny ]=size(A.FV.u);

    nt=length(got);
 got0 =  find(mod(mt(:),3/24)==0);
      kt=0;
 for jt=1:nt,...
     it=got(jt);
           kt=kt+1;
     FV.mt(kt,1)=mt(jt);
 for iv=1: length(vars);
    clear XX YY IN dim cv
    cv=char(vars(iv));
    if ~isfield(A.FV,(cv)), continue ; end
    [dim]= size(A.FV.(cv));
    dims(iv,1:length(dim))=dim;
      clear Z ZZ Vq
   if length(dim)<=2,
   FV.(cv)(kt,:)=A.FV.(cv)(it,:);
   end  
   if length(dim)==3,
   FV.(cv)(kt,:,:)=A.FV.(cv)(it,:,:);
   end  
    if length(dim)>=4,
      for iz=1:nz,
      FV.(cv)(kt,iz,:,:)=A.FV.(cv)(it,iz,:,:);
      end 
    end
  end; %iv
  
    if (mod(FV.mt(kt),3/24)==0) || jt==nt,...
    db=datestr(FV.mt(1),'yyyymmdd-HHMM');
    de=datestr(FV.mt(end),'yyyymmdd-HHMM');
    FVgrd_nm=['FV_Grd' db '_' de];
    save([path_matg  FVgrd_nm '.mat'],'FV','-v7.3');
    disp(['saved :'  FVgrd_nm ', nt=', num2str(length(FV.mt)),'; ' , datestr(now,31)]);
         kt=0;
         k=k+1;
      if jt<nt,...
       FV=[];     FV.mask=G.mask; FV.k=k; 
       FV.x=G.xq; FV.y=G.yq; FV.h=G.Vq; FV.siglev=G.siglev; 
      end
    end
 end;%jt 

end
toc
disp([ ' End at: '  datestr(now,31),' dt=' num2str((now-now1)*24*3600) ' s']); 
return
%% get the sigma-layers depth zs
    nt=length(FV.mt);
   zs=zeros(nt,nz+1,nx,ny);
for it=1:nt,...
 for i=1:nx,
 for j=1:ny,...
      yz =-(-(FV.h(i,j)+FV.el(it,i,j) ).*FV.siglev(1:end) ) + FV.el(it,i,j);
 zs(it,i,j,:)=yz;
 end; 
 end;
end

 