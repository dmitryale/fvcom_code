%function[]=plt_FV3_subset();
% call from Plot_FV_subset.m L.90
%(c) dmitry.aleynik@sams.ac.uk 2017.09.01                                 &
%  __o_O__¬                                                    FASTNEt    &
%  /_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS Glider mission 4, Malin Shelf &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
load(['../mat/' 'FV_3Layers_Grd_20170906-0810-20170907-1600.mat'],'-mat');
FV=FV3;

  this_dir=pwd; DRV=this_dir(1:1);
  path_fig='../fig_plt3/';
  path_mat='../mat/';
  [status,message,messageid]=mkdir(path_fig);   clear status message messageid
      Tb=datenum('2017-09-06 08:00:00',31); % dye release 1kg rodamin
     Te=datenum('2017-09-07 16:00:00',31);
addpath M:/Mar_Phys/matlab/m_map/
addpath M:/Mar_Phys/matlab/seawater/
addpath M:/Mar_Phys/matlab/general/

load('A_Grid.mat','-mat'); xv=G.xv;  yv=G.yv;  % x,y limits

if 0, % elevation timeseries
 filePO=['Poltips_Oban_20170901-0930.txt'];
 [PO] = import_pol_oban(filePO);
         cvar='el'; cunit='m';
    ist=10;    jst=10; % near south-westren border
    O.lon=FV.x(ist,jst);    O.lat=FV.y(ist,jst);
    O.ist=ist;  O.jst=jst; 
    
figure(1);clf;
cH2=plot(  FV.mt      ,FV.el(:,ist,jst),'-b.'); hold on;
cH1=plot(  PO.mtime   , PO.el,'k'); hold on;
 cH2.LineWidth=1.5;
 
% xlmt=[fix(min(FV.mt))-1 ceil(max(FV.mt))+0.1];
  xlmt=[fix(min(PO.mtime))+5.0 ceil(max(PO.mtime))+0.1];
  xlmt(2)=xlmt(1)+3.05;
set(gca,'xlim',xlmt);
cdates=[datestr(xlmt(1),'yyyymmdd') '-' datestr(xlmt(2),'yymmdd') ];

datetick('x','keeplimits')
  set(gca,'YMinorTick','on','XMinorTick','on');
 grid on
 hlg=legend([cH1,cH2],{['Poltips'];['FVcom'];});     
 set(hlg,'location','best','box','off');
title([ 'FV-grid ' cvar ', ' cunit ' at: ' num2str(O.lon) '^oW, ' num2str(O.lat) '^oN' ]);
         set(gcf,'PaperPositionMode','auto');

saveas(gcf,[path_fig 'FV_PO_el_dates' cdates ],'png')
end

if 1,...% bathymetry
Vq=[]; Vq=-FV.h; cvar='h'; cunit='m'; caxh=[-52 0];
it=1; Hmx=max(max(-Vq(isfinite(Vq))));

  figure(3); clf;
    hlv=[-250:25:-100 -90:10:-20 -15:5:-5 3 1 0 ];
    hmp=contourf(FV.x,FV.y,Vq,hlv); hold on;         
   caxis(caxh);
        plot(G.ncst(:,1),G.ncst(:,2),'.k-','markersize',0.5);
           daspect =[1  cosd(mean(yv)) 1];
  set(gca,'dataaspectratio',daspect);
  
  set(gca,'YMinorTick','on','XMinorTick','on');
  box on;
  set(gca,'xlim',xv);
  set(gca,'ylim',yv);
   ylabel('^oN');
   xlabel('^oW');

    cbar(1)=colorbar('location','east');
    set(cbar(1),'position',[ 0.84 0.20 0.013 0.65]);  % [ 0.053 0.18 0.02  0.56]);
    set(get(cbar(1),'title'),'string',[ cvar ',' cunit ]);
    set(gcf,'PaperPositionMode','auto');

   cb=datestr(FV.mt(it),'yyyymmdd-HHMM');
%  ce=datestr(Te,'yyyymmdd');
title([ 'FV-grid ' cvar ', ' cunit ', ' cb ]);
    set(gcf,'PaperPositionMode','auto');

           fig_names2  =([path_fig 'FV_grd_' cvar '_' cb  ]); % '_'  de ]);srid
print(['-f'],'-dpng','-loose','-r300',[fig_names2  'g.png']);
end
%% =========
%  vertical transect along Latitude or longitude = i
if 0,
[nt,nz,nx,ny] = size(FV.u);
%      iq=fix(nx/2); % half 
       iq=64; % Connel bridge
       O.Tlon=FV.x(iq,:);   O.Tlat=FV.y(iq,:);
                         clat=num2str(O.Tlat(1));
        it=1;
%       it=length(FV.mt);

 uu=squeeze(FV.u(it,:,iq,:)); 
 vv=squeeze(FV.v(it,:,iq,:));
 uv=sqrt(uu.^2+vv.^2); 
   Vq=uv; cvar ='uv'; cunit='m\cdots^{-1}'; 
           uvmx=max(max(abs(Vq(isfinite(Vq))))); caxf=[0.0 uvmx ];
%  Vq=u; cvar ='u'; cunit='m\cdots^{-1}'; caxf=[-0.3 0.3]; 
%  Vq=v; cvar ='v'; cunit='m\cdots^{-1}'; caxf=[-0.3 0.3];  

 el=squeeze(FV.el(it,iq,:));
 hh=squeeze(FV.h( iq,:));
 yz=zeros(size(uu));
%  to fix surface and seabed:
 for ii=1:length(hh)
  yz =-( -( hh(ii)   + el(ii )  ) .* FV.siglev(1:end-1) ) + el(ii);
% yz =-( -( FV.h(i,j)+FV.el(it,i,j) ).*FV.siglev(1:end) ) + FV.el(it,i,j);
 end

%yz=squeeze(FV.zs(iq,:,:))';
 xz=squeeze(FV.x(iq,:)); xz=repmat(xz,nz,1);
  
      figure(4);clf
  contourf(xz,yz,Vq); hold on
  pcolor(xz,yz,Vq); hold on; shading interp ; 
caxis([caxf]);
                    clv  =[-4:0.1:4]; 
                    clv05=[-4:0.5:4];
 [ hc, hp]=contour(xz,yz,Vq,[clv05]); hold on
    hp.LineColor='k';
    plot(xz,-hh,'-k'); hold on; % seabed line
    
    plot(xz,yz,'.k','markersize',1); hold on; % seabed line
    
 ylabel('Depth,m');
 xlabel('Longitude,^oW');
  db=datestr(FV.mt(it),'yyyymmdd-HHMM');
 set(gca,'YMinorTick','on','XMinorTick','on');

    cbar(1)=colorbar('location','east');
set(cbar(1),'position',[ 0.91 0.20 0.013 0.65]);  % [ 0.053 0.18 0.02  0.56]);
set(get(cbar(1),'title'),'string',[' ' cvar ',' cunit]);
 title([ 'FVcom: transect along ' clat '^oN ' db ', ' cvar ' '  cunit]);
  fig_name4  =([path_fig 'FV_grd_' cvar '_trans_' clat  'N_' db  ]);
  print(['-f'],'-dpng','-loose','-r300',[fig_name4  '.png']);
end;
%%

%% =========
 %% horizontal map:
if 1,... 
 figure(5); clf
  iz=2; it=1;
    Uq=squeeze(FV.u(it,iz,:,:));  caxf=[-4.0 4.0]; cvar ='u'; cunit='m\cdots^{-1}';
    Vq=squeeze(FV.v(it,iz,:,:));  caxf=[-4.0 4.0]; cvar ='v'; cunit='m\cdots^{-1}';
    UVq=sqrt(Uq.^2+Vq.^2);        caxf=[0.0  3.0]; cvar ='uv'; cunit='m\cdots^{-1}';
     Umx=max(max(abs(UVq(isfinite(UVq)))));
   
  pcolor(FV.x,FV.y,UVq); hold on; shading interp ; 
  caxis([caxf]);
    cbar(1)=colorbar('location','east');
set(cbar(1),'position',[ 0.92 0.20 0.013 0.65]);  % [ 0.053 0.18 0.02  0.56]);
set(get(cbar(1),'title'),'string',[' ' cvar ',' cunit]);
% colormap jet;
                            clv  =[-4:0.1:4];
                            clv05=[-4:0.5:4];
%      hmpf  = contourf(FV.x,FV.y,UVq,[clv  ]     ); hold on;
 [hmc, hmp ] = contour( FV.x,FV.y,UVq,[clv05],'-k'); hold on;
% clabel(hmc,[-4:1:4],'LabelSpacing',250);

plot(G.ncst(:,1),G.ncst(:,2),'.k-','markersize',0.5); hold on;
   set(gca,'ylim',[yv(1)+0.00015 yv(2)-0.00015]);
   set(gca,'xlim',[xv(1)+0.00015 xv(2)-0.00015]);
  box on;
title([ 'FV-grid ' cvar ', ' cunit ', iz=' num2str(iz) ', ' cb]);
   ylabel('^oN');
   xlabel('^oW');

   set(gca,'TickDir','out'); %in
   set(gca,'YMinorTick','on','XMinorTick','on');
   daspect4 =[1  cosd(mean(yv)) 1];
   set(gca,'dataaspectratio',daspect4);
try  plot(O.lon,O.lat,'rp','markersize',5);hold on; end;% elevation site
try  plot(O.Tlon,O.Tlat,'r-.'); hold on; end;          % transect line
       set(gcf,'PaperPositionMode','auto');

           fig_names5  =([path_fig 'FV_grd_' cvar 'iz' num2str(iz) '_' cb  ]);
 print(['-f'],'-dpng','-loose','-r300',[fig_names5  '.png']);

 end
