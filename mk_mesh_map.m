%% mk_mesh_map.m ; % (c)dmitry.aleynik@sams.ac.uk, 2013.05.30
% how to embedded  my AZIMUTH-FVCOM  mesh / data output into British coastline:
% - with m_map  
% - with the uk_coastline file 
%% ~~~ <^))))>< ~~~ ><(((@> ~~~

 clear all
 close all
% addpath(genpath([ 'M:\Mar_Phys\matlab\']));
  addpath( ([ 'M:\Mar_Phys\matlab\general']));
  addpath( ([ 'M:\Mar_Phys\matlab\images'])); % to access tricontour
  addpath( ([ 'M:\Mar_Phys\matlab\nan']));
 
 path_T1=['c:\HECTOR\matlab\azimuth3_results\'];  
  
        caseN ='3';  casename=['azimuth' caseN]; crun='v10_tmd\'; %crun='v9d_mf\';
       pathout=[casename '_' crun ];             
      FVCOM_plot_dir=[ path_T1 '\plots\' ];
      
    load([ FVCOM_plot_dir pathout 'mesh']);
      dir_lcvm = ls ([FVCOM_plot_dir pathout 'CVM_' '*.mat']);
                                  fmatc=dir_lcvm(1,:);    
     load([FVCOM_plot_dir pathout fmatc ]); 
      tri = mesh.trinodes(1,1:3); 
     tri_id=1; 
    xtri = mesh.nodexy(tri(1:3),1);
    ytri = mesh.nodexy(tri(1:3),2); 
    
     M=mesh.Nverts;
    N=mesh.Nelems;
 load('uk_coastline.mat');
l=find(isnan(coast(:,1)));

%% ==========
for  icon=1:1,...
%for  icon=2:4,...
       clear c h Clev  vrm varm
 %% a - mesh within m_map       
addpath( ([ 'M:\Mar_Phys\matlab\m_map']));
figure
  clf
  set(gcf,'position',[155  75 700   850]) ;        %  [55 55 1200   650]);  
  [IM_range_lat ,IM_range_lon] = deal([55.13  56.90], [ -7.10    -4.90]);   cm='';
  m_proj('UTM','lon',[ IM_range_lon],'lat',[IM_range_lat],'zon',30,'ell','wgs84') 
  m_grid('box','fancy')
  col=0.7*[1 1 1];
   [X,Y]=m_ll2xy(mesh.geog(:,1),mesh.geog(:,2),'clip','on');
   jnk=patch('Vertices',[X,Y],'Faces',mesh.trinodes,...
	   'EdgeColor',col,'FaceColor','none'); 
  view(-3,90)
  hold on
%% b - mesh without m_map       
figure
  clf
  set(gcf,'position',[55  55 700   850]) ; %  [55 55 1200   650]);

%% get the contour \\\\\\\\\\\\\\\\\\\\\\
 vrm=mesh.depth;
 vrm_i01=find(vrm<=0.01);
 mnf =nanmin( nanmin(vrm)' );
 mxf =nanmax( nanmax(vrm)' );  men =nanmean(nanmean(vrm)');
 mvs  =std(vrm);
 vrm(vrm<=0.01)=vrm(vrm_i01+1);
 
%    Nlev=15; %
     Clev1=[0:50:ceil(mxf)]; 
     [c1,h1] = tricontour(mesh.geog,mesh.trinodes,vrm,Clev1); hold on
    set(h1,'LineWidth',1.6,'linestyle','-');%
%   set(h1,'edgecolor','k','LineWidth',1.6,'linestyle','-');%
%   text_handles1 =  clabel(c1,h1,'LabelSpacing',125,'color','k','FontSize',12,'FontWeight','bold');
%% or  patch \\\\\\\\\\\\\\\\\\    
%   HPatch= patch('Vertices',[mesh.geog(:,1),mesh.geog(:,2)],'Faces',mesh.trinodes,...
%           'Cdata',vrm ,'edgecolor','interp','facecolor','interp');
        cbar=colorbar;  cbpos=[ 0.79    0.12    0.022    0.24];
    set(cbar,'position',cbpos);
%% or mesh     
%% \\\\\\\\\\\\\\\\\\
    
xlim2=[min(mesh.geog(:,1)) max(mesh.geog(:,1)) ];
ylim2=[min(mesh.geog(:,2)) max(mesh.geog(:,2)) ];
set(gca,'xlim',xlim2,'ylim',ylim2)
  
%                                 rlt=56;
                                  rlt=mean(mesh.geog(:,2));
set(gca,'dataaspectratio',[1 cosd(rlt) 1]); 
  hold on;
 
  ax1 = gca; 
  set(ax1,'box','on')
   
 xlim= get(gca,'xlim'); dx=(xlim(2)-xlim(1))*0.1; xt=xlim(1)+0.15*dx;
 ylim= get(gca,'ylim'); dy=(ylim(2)-ylim(1))*0.1; yt=ylim(2)-0.1*dy; %yt=ylim(1)+0.35*dy;
  
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4','PaperPosition',[.5 3 21 24]);
set(gcf,'PaperPositionMode','auto')
%   text_handles2 =  clabel(c2,h2,'LabelSpacing',145,'color','k','FontAngle','italic');% 144); % 72);
    hold on
 fc = get(gcf, 'Children');
  ylabel( ['{^o}N' ] );
 xlabel( ['{^o}W' ] );
 hold on;
 
  for i=1:size(l)-1
 	    fill(coast(l(i)+1:l(i+1)-1,1),coast(l(i)+1:l(i+1)-1,2),[.95 .95 .95])
  end
  
     title('AZIMUTH-FVCOM - bathymetry');
     
 set(gca,'xlim',xlim,'ylim',ylim)
 set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1])
 set(gcf,'PaperPositionMode','auto')
 hold on;
 
     fig_name=['FVCOM_mesh_map' ];
    namepng=[  fig_name '-bw'];
%   saveas(gcf,namepng,'fig')
   print(['-f'],'-dpng','-loose','-r300',[namepng '.png']);
   

end
%%