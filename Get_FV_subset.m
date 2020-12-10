%function[]=Get_FV_subset();
%(c) dmitry.aleynik@sams.ac.uk 2017.09.01                                 & 
%  __o_O__¬                                                    FASTNEt    &
%  /_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS Glider mission 4, Malin Shelf &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 clear all 
 close all
  this_dir=pwd; DRV=this_dir(1:1);
  path_fig='../fig/';
  path_mato='../mat/';
     Tb=datenum('2017-09-06 08:00:00',31); % dye release 1kg rodamin 
     Te=datenum('2017-09-07 16:00:00',31);
tic;
addpath M:/Mar_Phys/matlab/m_map
addpath M:/Mar_Phys/matlab/seawater
addpath M:/Mar_Phys/matlab/general/
  a=load('etive6_0coast.mat');

% path_mat=([DRV ':\HECTOR\matlab\etive8_results\plots\etive27_2014_e\']);
% load([path_mat 'FV_TS_20140205-20140208.mat']);
% path_mat=([DRV ':\HECTOR\matlab\etive8_results\plots\etive27_2013\']);
path_mat=([DRV ':\HECTOR\matlab\etive8_results\plots\etive27_2013_e\']);
         dir_fvts=dir ([path_mat 'FV_TS_20*.mat']);
 load([path_mat dir_fvts(end).name ]);
 FV_TS.names_station =FV_TS.name_station';

 filePO=['Poltips_Oban_20170901-0930.txt'];
 [PO] = import_pol_oban(filePO);
        mt_dif=PO.mtime(1)-FV_TS.mtime(1);       
        v_po = datevec(PO.mtime(1)   );
        v_fv = datevec(FV_TS.mtime(1));
        mt_difY=v_po(1)-v_fv(1); % 4 years 
        mt_difD=datenum(num2str(v_po(1)),'yyyy')-datenum(num2str(v_fv(1)),'yyyy'); %days
% Foram is ist=8 ; 
                          ist=8;
%                      mt_shft=mt_dif+5+9/24; %FV_TS_20140205;
%                      mt_shft=mt_difD-0; % 5+2/24; %FV_TS_20130620;
%                      mt_shft=mt_difD-0.75/24; 
                       mt_shft=mt_difD-0.50/24;  
if 1,
cH1=plot(  FV_TS.mtime+mt_shft,FV_TS.zeta(ist,:),'b'); hold on;
cH2=plot(  PO.mtime   , PO.el,'k'); hold on;
cH3=plot(  FV_TS.mtime+mt_shft,FV_TS.zeta(29,:),'g'); hold on;%ganavan

  xlmt=[fix(min(FV_TS.mtime+mt_shft))+0.0 ceil(max(PO.mtime))+0.1];
  xlmt=[fix(min(PO.mtime))+5.0 ceil(max(PO.mtime))+0.1];
 xlmt(2)=xlmt(1)+3.05;
set(gca,'xlim',xlmt);
cdates=[datestr(xlmt(1),'yyyymmdd') '-' datestr(xlmt(2),'yymmdd') ];

datetick('x','keeplimits')        
  set(gca,'YMinorTick','on','XMinorTick','on');
 grid on
     hlg=legend([cH1,cH2],{['Poltips'];['FVcom'];});
     set(hlg,'location','best','box','off');
 
  saveas(gcf,[path_fig 'PO_FV_el_dates' cdates ],'png')
end
  dir_fvm=dir ([path_mat 'FVCOM1_all*.mat']);
  
   load([path_mat dir_fvm(1).name ]);
%  load([path_mat 'FVCOM1_all140206_0206_11.mat']);
  try load([path_mat  'mesh.mat'],'-mat'); end;
  try load([path_mat 'CVM_090101.mat'],'-mat'); end
 time_offset = 678942; % from FVCOM time to matlab time
% [IM_range_lat ,IM_range_lon] = deal([56.433 56.50 ],[ -(5+28.0/60) -(5+22.5/60)]); cm='_AM'; %Ardmucknish
  [IM_range_lat ,IM_range_lon] = deal([56.445 56.496667 ],[ -(5+28.0/60) -(5+22.5/60)]); cm='_AM'; %Ardmucknish
   ml2m=1852*60; rlat= mean(IM_range_lat);
  xv=[IM_range_lon];   yv=[IM_range_lat];

%m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','grs80')
 m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','wgs84') 
 [mesh.XC,mesh.YC]=m_xy2ll(mesh.uvnode(:,1),mesh.uvnode(:,2));
if 1,...
 figure(2);clf
m_grid('box','fancy')
m_usercoast([ 'etive6_0coast.mat'],'Color','k','LineWidth',2);
%plot vertices
   col=0.7*[1 1 1];
   [X,Y]=m_ll2xy(mesh.geog(:,1),mesh.geog(:,2),'clip','on');
   jnk=patch('Vertices',[X,Y],'Faces',mesh.trinodes,...
	   'EdgeColor',col,'FaceColor','none');
view(-3,90);

nz=length(mesh.sigvec);
    %%!da get all cell, attached to that nodes
    
       [xx,yy,utmz] = deg2utm(56.5,-5.35);  %   utmz=['30 V'];  
       clear xx yy       
     if 0           
         nc=length(mesh.uvnode(:,1));
            utmzone(nc,1:4)=utmz;
       for  ic=1:nc, ...
            utmzone(ic,1:4)=utmz;
       end    
  [mesh.XC(:,1),mesh.YC(:,1)] = utm2deg( mesh.uvnode(:,1),mesh.uvnode(:,2),utmzone);  
     end
    
    X=mesh.geog(:,1); Y=mesh.geog(:,2); Z=mesh.depth;
  xv=[IM_range_lon]; yv=[IM_range_lat];
 [IN, ON]   = inpolygon(X,Y,xv,yv);
 [INC, ONC] = inpolygon(mesh.XC,mesh.YC,xv,yv);
%  ntr=60; %7.4km x 5.7km
%  dvx=(xv(2)-xv(1))/ntr;  %in degrees
%  dvy=(yv(2)-yv(1))/ntr;
      dv=20;%m
  dvy=dv/ ml2m ; 
  dvx=dv/(ml2m * cosd(rlat));
 xx=[xv(1):dvx:xv(2)]';  yy=[yv(1):dvy:yv(2)]';    
     XX=X(IN);    YY=Y(IN);    ZZ= Z(IN);
     XXc=mesh.XC(INC);  YYc=mesh.YC(INC);  
        [xq,yq] = meshgrid(xx,yy );     % plot(XX,YY,'.'); plot(xq,yq,'.'); hold on 
[Xq,Yq,Vq]  = griddata(XX,YY,ZZ,xq,yq) ; % , method) triscateredinterp
 [INb, ONb] = inpolygon(Xq,Yq,a.ncst(:,1),a.ncst(:,2));
Vq0=Vq;
Vq(~INb)=NaN;
mask=zeros(size(Vq)); 
mask(Vq>0)=1; maskn=mask; maskn(mask<=0)=NaN;
%% ======================================
cels=[];
in=find(IN>0);
for ii=1:length(in),
    no=in(ii);
    ce=[];ceu=[];
    for jj=1:3,
        cel=find(mesh.trinodes(:,jj)==no);
    ce(jj,1:length(cel))=cel;
    end  
                    ce=ce(ce>0);
         ceu=unique(ce);         
    cels=cat(1,cels,ceu);
    
    nc=length(ceu);
    G.no(ii,1)=no;
    G.nc(ii,1)=nc;
    G.ceu(ii,1:nc)=ceu;
      lonc=[]; latc=[]; zu=[];
    for jj=1:nc,...
         cenj=mesh.trinodes(ceu(jj),:);    
     lonc(jj,1)=mean(mesh.geog(cenj,1));
     latc(jj,1)=mean(mesh.geog(cenj,2));
     for ik=1:nz,...
       zu(jj,ik)=mesh.zuv(ik, ceu(jj) );
     end
    end
    G.lonc(ii,1)=mean(lonc);
    G.latc(ii,1)=mean(latc);
    G.depc(ii,1)=mean(mesh.uvdepth(ceu));
     
for ik=1:nz,...
    G.zu(  ii,ik)=mean(zu(:,ik));
end 

    G.lon(ii,1) = mesh.geog(no,1);
    G.lat(ii,1) = mesh.geog(no,2);
    G.dep(ii,1) = mesh.depth(no);
    G.z(  ii,:)  = mesh.nodez(:,no);
    
end
cels=unique(cels);        
G.mask=mask;
G.maskn=maskn;
G.inc=cels; 
G.xq=xq;
G.yq=yq;

G.Xq=Xq;
G.Yq=Yq;
G.Vq=Vq;
G.IN =IN;
G.INC=INC;
G.INb=INb; % outside boundaries
 G.XX= XX;   G.YY= YY ;
 G.XXc= XXc;  G.YYc= YYc ;
G.xv=xv; G.yv=yv;
G.xx=xx; G.yy=yy;
G.utmz=utmz;
G.ncst=a.ncst;
G.dv=dv; G.dvx=dvx; G.dvy=dvy;
G.siglay=double(FV_TS.siglay(1,:)');
G.siglev=double(FV_TS.siglev(1,:)');
save(['A_Grid.mat'],'G','-v7.3');

figure(2);clf
        plot(G.xq(G.INb),G.yq(G.INb),'.');hold on;
        plot(G.ncst(:,1),G.ncst(:,2),'.k-','markersize',0.5);        
   set(gca,'YMinorTick','on','XMinorTick','on');             
           daspect =[1  cosd(mean(yv)) 1];
  set(gca,'dataaspectratio',daspect);      box on;
  set(gca,'xlim',xv);  set(gca,'ylim',yv);

else
    load('A_Grid.mat','-mat');
end
%% ======================

if 0,...
   got= find(FVCOM1.Date+mt_shft >=Tb & FVCOM1.Date+mt_shft <=Te);
it=got(1);
     Z = FVCOM1.el(it,:); cvar ='el'; cunit='m'; % 'cm^2s^{-2}';
     ZZ = double(Z(G.IN));
  [Xq,Yq,Vq] = griddata(G.XX,G.YY,ZZ,G.xq,G.yq) ; % , method) triscateredinterp
  Vq0=Vq; Vq(~G.INb)=NaN;  
  
  figure(3); clf;
  contourf(Xq,Yq,Vq); hold on;
    plot(G.ncst(:,1),G.ncst(:,2),'.k-','markersize',0.5);
%   plot(G.lon,G.lat,'.k'); hold on;
%   plot(G.lonc,G.latc,'.r'); hold on

   set(gca,'YMinorTick','on','XMinorTick','on');
             
           daspect =[1  cosd(mean(yv)) 1];
  set(gca,'dataaspectratio',daspect);    
  box on;
     
  set(gca,'xlim',xv);
  set(gca,'ylim',yv);
    cbar(1)=colorbar('location','east'); 
    set(cbar(1),'position',[ 0.84 0.20 0.013 0.65]);  % [ 0.053 0.18 0.02  0.56]);  
    set(get(cbar(1),'title'),'string',[ cvar 'w,' cunit ]);
    set(gcf,'PaperPositionMode','auto');
  
   cb=datestr(FVCOM1.Date(it)+mt_shft,'yyyymmdd-HHMM');
%  ce=datestr(Te,'yyyymmdd');
title([ 'FV-grid ' cvar ', ' cunit ', ' cb ]);

           fig_names2  =([path_fig 'FV_grd_' cvar '_' cb  ]); % '_'  de ]);srid
print(['-f'],'-dpng','-loose','-r300',[fig_names2  'g.png']);
end
%% =========
now1=now;
  disp([ 'start at: '  datestr(now1,31)]);
 [fn fi]=size(dir_fvm);
    k=0;
for fi=1:fn,...
   load([path_mat dir_fvm(fi).name ]);
   got=[];
   got= find(FVCOM1.Date+mt_shft >=Tb & FVCOM1.Date+mt_shft <=Te);
   if isempty(got), continue;end
 
%% define vertical levels:
   FV=[]; FV.mask=G.mask;
FV.x=G.xq; FV.y=G.yq; FV.h=G.Vq;   FV.siglev=G.siglev;
[nx,ny  ]=size(FV.h);
[kt,m]   =size(FVCOM1.el);
[kt,nz,m]=size(FVCOM1.t);
[kt,nz,n]=size(FVCOM1.u);
% nz=length(mesh.sigvec);

             nt=length(got);
             mt0=FVCOM1.Date(got);
             mt=FVCOM1.Date(got)+mt_shft;
       got0=  find(mod(mt(:),3/24)==0);
       
            cvar ='el'; cunit='m'; % 'cm^2s^{-2}';
   vars={'el';'t';'s'; 'u';'v';'ww';'ua';'va';'uwind_speed';'vwind_speed'};
%  vars={'el';'t';'u';};
   FV.x=G.xq; FV.y=G.yq; FV.h=G.Vq; FV.siglev=G.siglev; %FV.zs=zs;
%  
     kt=0;
 for jt=1:nt,...
     it=got(jt);
           kt=kt+1;
     FV.mt(kt,1)=mt(jt);
 for iv=1: length(vars);
  clear XX YY IN dim cv
    cv=char(vars(iv));
    [dim]= size(FVCOM1.(cv));
    dims(iv,1:length(dim))=dim;
  if dim(end)==mesh.Nverts,...
      XX=G.XX;  YY=G.YY;   IN=G.IN;
  else
      XX=G.XXc;  YY=G.YYc; IN=G.INC;
  end
      clear Z ZZ Vq
   if length(dim)<=2,
       Z = FVCOM1.(cv)(it,:);
      ZZ = double(Z(IN));     
   [Xq,Yq,Vq] = griddata(XX,YY, ZZ, G.xq, G.yq) ; % , method) triscateredinterp
                  Vq(~G.INb)=NaN;
   FV.(cv)(kt,:,:)=single(Vq);
   end
    if length(dim)>=3,
      Z =  double(squeeze( FVCOM1.(cv)(it,:,:)));
     for iz=1:nz,
     ZZ = (Z(iz,IN))';     
     [Xq,Yq,Vq] = griddata(XX,YY, ZZ, G.xq, G.yq) ; % , method) triscateredinterp
        Vq0=Vq;          Vq(~G.INb)=NaN;
     FV.(cv)(kt,iz,:,:)=single(Vq);
     end 
    end
  end; %iv
  
 if 0,...
  for i=1:nx,
   for j=1:ny,...
     yz =-( -( FV.h(i,j)+FV.el(kt,i,j) ).*FV.siglev(1:end) ) + FV.el(kt,i,j);
   FV.zs(kt,:,i,j)=yz;
   end; 
  end;
 end
 
    if (mod(FV.mt(kt),3/24)==0) || jt==nt,...
    db=datestr(FV.mt(1),'yyyymmdd-HHMM');
    de=datestr(FV.mt(end),'yyyymmdd-HHMM');
    FVgrd_nm=['FV_Grd' db '_' de];
    save([path_mato FVgrd_nm '.mat'],'FV','-v7.3');
    disp(['saved :'  FVgrd_nm ', nt=', num2str(nt),'; ' , datestr(now,31)]);
         kt=0;
         k=k+1;
      if jt<nt,...
       FV=[];  FV.mask=G.mask; FV.k=k; 
       FV.x=G.xq; FV.y=G.yq; FV.h=G.Vq; FV.siglev=G.siglev;  
      end
    end
 end;%jt 

end
toc
disp([ ' End at: '  datestr(now,31),' dt=' num2str((now-now1)*24) ' h']); 

 if 1,
     iq=fix(nx/2);
     %  vertical transect along Latitude or longitude = i
       a=  FV.y(iq,:); clat=num2str(a(1));
                 it=length(FV.mt);
 uu=squeeze(FV.u(it,:,iq,:)); caxf=[-0.3 0.3]; cvar ='u'; cunit='m\cdots^{-1}';   
 el=squeeze(FV.el(jt,iq,:));
 hh=squeeze(FV.h(iq,:));
 yz=zeros(size(uu));
%  to fix surface and seabed:
 for ii=1:length(hh)
   yz(:,ii) =-( -( hh(ii)   + el(ii )  ) .* FV.siglev(1:end-1) ) + el(ii);
 end
   
%yz=squeeze(FV.zs(iq,:,:))';
 xz=squeeze(FV.x(iq,:)); xz=repmat(xz,nz,1);
      figure(4);clf
contourf(xz,yz,uu); hold on
caxis([caxf]);
           clv=[-1:0.05:1];
 [ hc hp]=contour(xz,yz,uu,[clv]); hold on
 hp.LineColor='k';
    plot(xz,-hh,'-k'); hold on
 ylabel('Depth,m');
 xlabel('Longitude,^oW');
  db=datestr(FV.mt(it),'yyyymmdd-HHMM');
 set(gca,'YMinorTick','on','XMinorTick','on');       
 
    cbar(1)=colorbar('location','east');
set(cbar(1),'position',[ 0.91 0.20 0.013 0.65]);  % [ 0.053 0.18 0.02  0.56]); 
set(get(cbar(1),'title'),'string',[' ' cvar ',' cunit]);
 title([ 'FVcom: transect along ' clat '^oN ' db ', ' cvar ' '  cunit]);
  fig_name4  =([path_fig cvar '_trans_' clat  'N_' db  ]);  
  print(['-f'],'-dpng','-loose','-r300',[fig_name4  '.png']);
 %%
 
 %% maps:
 
 figure(5); clf

     
 iz=1;
    Vq=squeeze(FV.u(it,iz,:,:));  caxf=[-4.0 4.0]; cvar ='u'; cunit='m\cdots^{-1}';   
 %  Vq(~G.INb)=NaN;  
 hmp=  contourf(Xq,Yq,Vq); hold on;
   caxis([caxf]);
           clv=[-4:0.1:4];
 
    cbar(1)=colorbar('location','east');
set(cbar(1),'position',[ 0.92 0.20 0.013 0.65]);  % [ 0.053 0.18 0.02  0.56]); 
set(get(cbar(1),'title'),'string',[' ' cvar ',' cunit]);
colormap jet
plot(G.ncst(:,1),G.ncst(:,2),'.k-','markersize',0.5); hold on;
   set(gca,'ylim',[ylm4(1)+0.0015 ylm4(2)-0.0015]);   
   set(gca,'xlim',[xlm4(1)+0.0015 xlm4(2)-0.0015]);

  box on;
title([ 'FV-grid ' cvar ', ' cunit ', ' cb ,]);
   do_gmap='Y';  % do_gmap='N';
if do_gmap=='Y',...
   ylm4=yv;  xlm4=xv;
  mtype={'roadmap','satelite','terrain', 'hybrid'};  gk=3;
  APIKey='AIzaSyC55hr1J9e8y1LpbjTSUaogudBKRc_X4uo';
  api_key=APIKey; apiKey=APIKey; 
          scale=2;       
          mtype_g = char(mtype(gk));
  clear  hg  %[lonVec3 latVec3 imag3] 
  plot(xv,yv,'w.'); hold on;
 hg= plot_google_map( ... % 'height',ny,'width',nx,...  'height',640,'width',640,...
 'maptype',mtype_g,... 
 'alpha',0.85 ,...  % (0-1) Transparency level of the map (0 is fully transp)
 'refresh',1,...
 'autoAxis',0,...
 'showLabels',0,...
 'language','',...
 'marker','yellowA');
   set(gca,'ylim',[ylm4(1)+0.0015 ylm4(2)-0.0015]);   
   set(gca,'xlim',[xlm4(1)+0.0015 xlm4(2)-0.0015]);
end
   set(gca,'TickDir','out'); %in
   set(gca,'YMinorTick','on','XMinorTick','on');
   daspect4 =[1  cosd(mean(ylm4)) 1];
   set(gca,'dataaspectratio',daspect4);
           fig_names5  =([path_fig 'FV_grd_' cvar '_' cb  ]);  
print(['-f'],'-dpng','-loose','-r300',[fig_names5  'g0.png']);

 end
 