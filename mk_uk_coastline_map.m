%% mk_ukcoastline_map.m ; % (c)dmitry.aleynik@sams.ac.uk, 2013.05.30
%  British coastline:
% - with m_map  
% - with the uk_coastline file 
%% ~~~ <^))))>< ~~~ ><(((@> ~~~

 clear all
 close all
  addpath( ([ 'M:\Mar_Phys\matlab\general']));
  addpath( ([ 'M:\Mar_Phys\matlab\nan']));
%  path_T1=['c:\HECTOR\matlab\azimuth3_results\'];  
load('uk_coastline.mat');
l=find(isnan(coast(:,1)));

%% ==========
       clear c h Clev  vrm varm
addpath( ([ 'M:\Mar_Phys\matlab\m_map']));
%% b - mesh without m_map
figure
  clf
  set(gcf,'position',[55  55 700   850]) ; %  [55 55 1200   650]);
%    Nlev=15; %
  mxf=350 ;% max depth
%      Clev1=[0:50:ceil(mxf)]; 
%         cbar=colorbar;  cbpos=[ 0.79    0.12    0.022    0.24];
%     set(cbar,'position',cbpos);
%% or mesh     
%% \\\\\\\\\\\\\\\\\\
ylim2=[ 55.13  56.90];
xlim2=[ -7.10  -4.90];
set(gca,'xlim',xlim2,'ylim',ylim2)
  
                                  rlt=56;
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
  
     title('uk - bathymetry');
     
 set(gca,'xlim',xlim,'ylim',ylim)
 set(gca,'dataaspectratio',[1 cosd(mean(ylim)) 1])
 set(gcf,'PaperPositionMode','auto')
 hold on;
 
     fig_name=['uk_coast_map' ];
    namepng=[  fig_name '-bw'];
   print(['-f'],'-dpng','-loose','-r300',[namepng '.png']);
   

%%