%function[]=plt_FV_subset_vect_g();
%(c) dmitry.aleynik@sams.ac.uk 2017.09.01                                 &
%  __o_O__¬                                                    FASTNEt    &
%  /_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS Glider mission 4, Malin Shelf &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  this_dir=pwd; DRV=this_dir(1:1);
  path_fig='../fig_plt/';
  path_mat='../mat/';
  [status,message,messageid]=mkdir(path_fig);   clear status message messageid
     Tb=datenum('2017-09-06 08:00:00',31); % dye release 1kg rodamin
     Te=datenum('2017-09-07 16:00:00',31);
addpath M:/Mar_Phys/matlab/m_map/
addpath M:/Mar_Phys/matlab/seawater/
addpath M:/Mar_Phys/matlab/general/

% load('A_Grid.mat','-mat'); xv=G.xv;  yv=G.yv;  % x,y limits

%% =========
 %% horizontal map:
  xvec=xv(2)-0.6/60; yvec=yv(2)-0.6/60;
  scalev = 0.5; % arrow scale
  scale  = 5.0;
  rlat=mean(yv);
        
if 1,... 
 figure(6); clf 
 set(gcf,'position',[50 50 580 580]);
 colormap jet;
 
 iz=1; it=1;  
clear x y xx yy ug vg Uq Vq UVq go u v uv vp goa ;
% skip ploting magnitude uv map (subset with sps )   
if 0,...
           caxf=[-4.0 4.0];  cunit='m\cdots^{-1}';
    Uq=squeeze(FV.u(it,iz,:,:)); cvar ='u';   caxf=[-4.0 4.0];
    Vq=squeeze(FV.v(it,iz,:,:)); cvar ='v';   caxf=[-4.0 4.0];
    UVq=sqrt(Uq.^2+Vq.^2);       cvar ='V';   caxf=[0.0  3.0];    
     Umx=max(max(abs(UVq(isfinite(UVq)))));
   
  pcolor(FV.x,FV.y,UVq); hold on; shading interp ; 
  caxis([caxf]);
    cbar2=colorbar('location','west');
set(cbar2,'position',[ 0.80 0.33 0.013 0.35]);  % [ 0.053 0.18 0.02  0.56]);
set(get(cbar2,'title'),'string',[' ' cvar ',' cunit]);
 colormap jet;
                            clv  =[-4:0.1:4];
                            clv05=[-4:0.5:4];
%       hmpf  = contourf(FV.x,FV.y,UVq,[clv  ]     ); hold on;
  [hmc, hmp ] = contour( FV.x,FV.y,UVq,[clv05],'-k'); hold on;
% clabel(hmc,[-4:1:4],'LabelSpacing',250);
end;
%% ==============================
clear  Uq Vq UVq ig jg go 
        sps=4; % sps=1; % skipping points 
  ig=[1:sps:nx];
  jg=[1:sps:ny];   
                  caxf=[-4.0 4.0]; cvar ='u'; cunit='m\cdots^{-1}';
    Uq=squeeze(FV.u(it,iz,ig,jg)); 
    Vq=squeeze(FV.v(it,iz,ig,jg)); cvar ='v';  
    UVq=sqrt(Uq.^2+Vq.^2);         cvar ='V';   caxf=[0.0  3.0];     
     Umx=max(max(abs(UVq(isfinite(UVq)))));
     masks = squeeze(FV.mask(ig,jg));
           x=squeeze(FV.x(ig,jg)); 
           y=squeeze(FV.y(ig,jg)); 
  go=[]; go=find(masks >0); 
     xx=squeeze(x(go)) ;
     yy=squeeze(y(go)) ;         
     u=squeeze(Uq(go)) ;
     v=squeeze(Vq(go)) ;
     uv=squeeze(UVq(go)) ;
     u(abs(u)>5)=NaN; v(abs(v)>5)=NaN;
                   clv=[-4 : 0.2 : 4   ]; 
     
     nk=length(xx);
  xx(nk+1)=xvec;  yy(nk+1)=yvec; 
  xx(nk+2)=xvec;  yy(nk+2)=yvec;
   u(nk+1)=scalev;   v(nk+1)=0;
   u(nk+2)=0;        v(nk+2)=scalev; 
                     vp=v*cosd(rlat); % adjustment for plotting vectors 
   set(gca,'xlim',[xv]);set(gca,'ylim',[yv]);
   quiverc(xx,yy,u,vp,scale); hold on;
   caxis(caxf);
%  set(gcf, 'InvertHardCopy', 'off')
 hts=text(xvec(1)-0.00045 , yvec(1)-0.00075, [num2str(scalev) ' ms^{-1}'],'FontSize',9); hold on;

%% create Vectors' colorbar & ticks , that is not very accurate - only for plotting 
   c=[0: 0.5: max(caxf)];
   cbar=colorbar('location','EastOutside','XTickLabel',c);
                       cbposi = [0.6895 0.048 0.28 0.015];
                       cbposi = [ 0.92 0.20 0.013 0.65]; 
   set(cbar,'position',cbposi);      
%     ct=cbar.Ticks; cbtl=cbar.TickLabels;
if  cbar.Limits(end)<10,
    cbar.Ticks=[0:0.5:cbar.Limits(end)];
else
    cbar.Ticks=[0:10:cbar.Limits(end)];    
end
      cbar.TickLabels=c;
%% -----------------------      
      
   colormap jet;
%  freezeColors;
   set(gca,'xlim',[xv]);set(gca,'ylim',[yv]);

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
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'color','w'); % 'none')
    orient portrait ;

                 fig_names6  =([path_fig 'FV_Vect_grd_' cvar 'iz' num2str(iz) '_' cb  ]);
 
print(['-f'],'-dpng','-loose','-r300',[fig_names6  '.png']);

 end
