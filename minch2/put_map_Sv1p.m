% function[]=put_sst_map_Sv1p
%(c)dmitry.aleynik@sams.ac.uk, 2013.05.15)
%  __o_O__¬                                                               &
%  \_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS ~~~~~~~~~~                    &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
if 1,...
  clear;
  a=pwd;  DRV=a(1:1); this_dir = a;
  if ispc, ... % isunix
%          mhome='M:/';  
           mhome='./';
   path_omat =[DRV ':/samhanach/sa01da/work/minch2/Archive/'];
  else
             mhome='/home/sa01da/';
  path_omat=[mhome '/work/minch2/Archive/'];
  end             
%  addpath([mhome 'Mar_Phys/matlab/t_tide/'  ]);
   addpath([mhome 'Mar_Phys/matlab/seawater/']);
   addpath([mhome 'Mar_Phys/matlab/general/' ]);
   addpath(['extLib/']);

 data_dir =[this_dir '/test/'];  
 path_fig =[data_dir '/Figures/Maps/'];
 path_out = path_fig;
  [status,message,messageid] = mkdir(path_fig);  
  
%   load([mhome 'Mar_Phys/matlab/m_map/uk_coastline.mat']);
    load([  'uk_coastline.mat']);
  l=find(isnan(coast(:,1)));
 
 
  case_1=1;          dir_fvm=dir([path_omat '/mat_20*']);
%                      ny=1; %2013_a
                       ny=5; %2016
   path_M=[path_omat dir_fvm(ny,1).name '/'];
                      dir_mat =dir([path_M '/FV*_all*.mat']);
    [ nf,  lf]= size(dir_mat);
   for lf=1:nf,
         fvnm=[]; fvnm=dir_mat(lf,1).name;
         dir_m(lf,1:length(fvnm))=fvnm;
      mt_f(lf,:)= dir_m(lf,11:11+5);   
   end
   mtf=datenum(mt_f,'yymmdd');
           cyear=datestr(mtf(1),'yyyy');    %       cyear='2016'; 
           iyear=str2double(cyear);
Eom=[];  Eom = cumsum(eomday(iyear,[1:12])); mjj=[1:12];

% mj=[24*(Eom(mjj(8)) + 18)+1 : 24*(Eom(mjj(8)) +18)+ 24*1 + 1]; ite=max(mj); %spring tide 22 sept 2016
%  mtd=datenum('2016','yyyy')+fix(mj(1)/24); datestr(mtd)
%  [ df, kf]= min(abs((mtf-mtd)));
%       load([path_M '/' 'Mesh.mat']);
      load(['Mesh.mat']); 
kf=1;
%        nfiles=length(dir_mat);
         nfiles=1;
for kf=1:nfiles,...
FVCOM1=[]; FV=[];
            fl = dir_mat(kf,1).name ; 
           load([ path_M fl ],'-mat');
%% do SV           
      [nt,nz,N]=size(FVCOM1.t); % N=Mesh.Nverts;
  h=double(FVCOM1.h);
  KBM1=nz;   dsig=1/(KBM1+1);     sigvec=(0:(KBM1-1))'*dsig+0.5*dsig;
  FV=[]; FV.mtime=FVCOM1.Date; FV.sigvec=sigvec; FV.el=FVCOM1.el;
    ZI=[0,10,20,30,40 ,75,100,150,200,250];NZ=length(ZI);
    
for it=1:nt,
    clear t s r z Sv zeta  Svmx Svmx_zi Svmx_z
       t=squeeze(double(FVCOM1.t(it,:,:)));
       s=squeeze(double(FVCOM1.s(it,:,:)));
       r=squeeze(double(FVCOM1.r(it,:,:)));      
       zeta=squeeze(double(FVCOM1.el(it,:)))' ; 
       Sv=zeros(size(s));     
        z=zeros(nz,N); p=z;
       for ik=1:nz,...
        z(ik,:) =-sigvec(ik)*(h+zeta);
        p(ik,:) = sw_pres(-squeeze(z(ik,:)),Mesh.geog(:,2)')';% .lat);
       end
           Sv=sw_svel( s, t, p);
       FV.t(it,:,:) =t; 
       FV.s(it,:,:) =s; 
       FV.r(it,:,:) =r; 
      FV.p(it,:,:) =p; 
      FV.z(it,:,:) =z; 
      FV.Sv(it,:,:)=Sv; 
 % find the maximum & its z-index & depth ;
   [Svmx Svmx_zi]=max(Sv); 
   for ii=1:N
    Svmx_z(ii)=-z(Svmx_zi(ii),ii);
   end
    FV.Svmx(it,:)=Svmx ;
    FV.Svmx_z(it,:)=Svmx_z ;
    FV.Svmx_zi(it,:)=Svmx_zi ;
 %% reinterpolat at fix ZIs
   SvIZ=zeros(NZ,N);
   for ii=1:N,...
    Svi=[];   Svi=interp1([0;-z(:,ii)],[Sv(1,ii); Sv(:,ii)], ZI);
   SvIZ(:,ii)=Svi;
   end
   FV.SvIZ(it,:,:)=SvIZ;
end
  FV.ZI=ZI;
end
%% =========
% IZ=1;
  IZ=3; %30m
% IZ=5; %40m
%       no=length(G.obc_nodes);ni=fix(no/2);
%   w_elm=G.obc_elems(ni); w_nod=G.obc_nodes(ni);

do_log='Y'; do_log='N';
mkgif='Y'; gi=0; gim=0; 
%mkgif='N';

clear kt tml tb te jw mt jjt;

pos(1,:)=[0.06 0.69 0.41 0.30]; 
pos(2,:)=[0.06 0.36 0.41 0.30];
pos(3,:)=[0.06 0.045 0.41 0.30];
pos(4,:)=[0.53 0.69 0.41 0.30];
pos(5,:)=[0.53 0.36 0.41 0.30];
pos(6,:)=[0.53 0.045 0.41 0.30];
poss= [0.06 0.045 0.90 0.90]; 
%   ExpLim=0.0005; cm=''; %0.0005=0.05%
%   ExpLim=0.0002; cm=''; %0.0002=0.02%
    ExpLim=0.001 ; cm=''; %0.0001=0.01%
%   vs5=[-0:0.5: 42]; vs2=[-0:0.25:42];
%     vd05=[ExpLim,10*ExpLim, 0.01            , 0.1:0.05:1 ];  
%     vd01=[ExpLim,10*ExpLim, 0.01: 0.01: 0.09, 0.1:0.05:1 ];  
%  caxD =[ExpLim 1];
   
   caxf0=[1496, 1508 ];    vd01=[ExpLim,caxf0(1): 1 :caxf0(2)  ]; 
%  caxf0=[1490, 1502 ];    vd01=[ExpLim,caxf0(1): 1 :caxf0(2)  ]; 
%   caxf0=[0,  255 ];    vd01=[       caxf0(1): 10 :caxf0(2) ]; 
  caxD= caxf0; 
if 1,...
%     scale = 0.25;     scal_Q= 10; 
   y1 = 56.100; y2 = 56.87;  x1 = -6.40; x2 = -5.00; cm='FoL'  ; % F of Lorn
GA.rlon =[x1 x2];  
GA.rlat =[y1 y2]; 
 y1=min(GA.rlat)-0.02; y2=max(GA.rlat)+0.02; x1=min(GA.rlon)-0.2; x2=max(GA.rlon)+0.2;  
  
   xlm4=[x1 x2]; ylm4=[y1 y2]; xlm=xlm4; ylm=ylm4;
[IM_range_lat ,IM_range_lon] = deal(ylm4,xlm4); cm='zoo'  ;
 end
% try load([data_dir  'CTDxy160728-160802.mat']);end
% IZ=1;
%   cax8=[11.0 17];   %25 temp    
%   caxS=[23.0 35.5]; %25 sal
%     caxS=[31.0 36]; %25 sal
  do_gmap='N'; %'Y';
  mtype={'roadmap','satelite','terrain', 'hybrid'};
  gk=3; %'terrain', 
% gk=4; %'hybrid' % as chla
 
%   y1 = 54.5; y2 = 56.1; x1 = -6.2; x2 = -4.4;
%  xlm4=[x1 x2]; ylm4=[y1 y2];  
   xlm=xlm4; ylm=ylm4;
 clear igd; 
  igd = find (  Mesh.geog(:,1)<IM_range_lon(2) & Mesh.geog(:,1) > IM_range_lon(1) &...
                Mesh.geog(:,2)<IM_range_lat(2) & Mesh.geog(:,2) > IM_range_lat(1));
  jp6=0;              
         ie=length(FVCOM1.Date); 
         ie=1;
         Date=FVCOM1.Date;
for it=1:ie,
 %        [jk, jjt]=min(abs(Date(it)-FVCOM1.Date));
   if ( mod(it-1,1 )~=0), continue; end
%  if ( mod(it-1,4 )~=0), continue; end
 disp(['S.V.: it=', num2str(it) ', IZ=' num2str(IZ) ', ' datestr((FVCOM1.Date(it)),' yyyymmdd-HH') ' jp6=' num2str(jp6) ]);
%             ciz=num2str(IZ); 
              ciz=[num2str(ZI(IZ)) 'm'];  
            % ciz='Zmx';% Svmx
               jjt=it;
% if jp6==0, 
    figure(3000); clf; clear sp                
 %  set(gcf,'position',[50 50 520 580]);   
 %  set(gcf,'position',[50 50 580 540]);   
    set(gcf,'position',[50 50 680 800]);   
%end
      jp6=jp6+1;
      pos6=[];  pos6=pos(jp6,:);        pos6=poss;
   %sp(jp6)=subplot('position',pos6); 
   
   clear isf varm varn varmp ax3 HPatch

%  varm(:,:)=squeeze(FVCOM1.DYE(jjt,IZ,: ));
% varm(:,:)=squeeze(FV.Sv(jjt,IZ,: ));
% varm(:,:)=squeeze(FV.Svmx(jjt,: ))';
% varm(:,:)=squeeze(FV.Svmx_z(jjt,: ))';
  varm(:,:)=squeeze(FV.SvIZ(jjt,IZ,: ));
 
           varm(varm < 0.01*ExpLim)=0;   varm(varm <= 0)=0; %ExpLim*1.0e-10; % 5*10-4   
   varmp=double(varm);
%  varn(:,:)=squeeze(FVCOM1.DYE(jjt,IZ,igd ));
%    varn(:,:)=squeeze(FV.Sv(jjt,IZ,igd ));
%    varn(:,:)=squeeze(varm(igd ));
   varn(:,:)=squeeze(  varm(igd ));
   
              varn(varn < 0.1*ExpLim)=NaN;         
 isf=find(isfinite(varn)>0); varn=varn(isf);
 mnf=nanmin( nanmin(varn)' );
 mxf=nanmax( nanmax(varn)' ); 
 mef=nanmean(nanmean(varn)');
 mvrf=var(varn);
 mvsf=std(varn);
      caxf1=[mnf mxf];
 caxD=caxf0;
         wst=varmp;
         lgvd01   =vd01;
  if do_log=='Y';         
          wst=log10(wst);          
       lgvd01=[10.^(log10(ExpLim):log10(-1))];
    wslog_mx=[min(min(min(wst))), max(max(max(wst)))];
  end;
  
   HPatch=patch('Vertices',Mesh.nodexy,'Faces',Mesh.trinodes,... 
      'Cdata',wst ,'edgecolor'   ,'interp','facecolor','interp'); hold on
% %    'Cdata',varmp ,'edgecolor','interp','facecolor','flat'); hold on
% % [c,h] =  
      tricontour(Mesh.geog,Mesh.trinodes,wst,lgvd01); hold on;
      
% M=length(Mesh.geog);
% H=tricontour_delf(Mesh.trinodes,Mesh.geog(:,1),Mesh.geog(:,2),zeros(M,1),wst,vd01); hold on;
% H=tricontour_delf(tri,x,y,z,v,levels,color)
%   
%   [CS,h]=tricontf(Mesh.geog(:,1),Mesh.geog(:,2),Mesh.trinodes,wst);  hold on
%   set(h,'edgecolor','none'); hold on;
%          tricont(Mesh.geog(:,1),Mesh.geog(:,2),Mesh.trinodes,wst,'-k'); hold on;
%
% Nlev=15; Clev=[0:0.1:3, 3:0.5:42];  
%[c,h] = tricontour([x(:,1),x(:,2)],tri,x(:,3),Clev);
%[c,h] = tricontour(Mesh.nodexy,Mesh.trinodes,varm,vs2 ); hold on;% Clev);
% clabel(c,[0:1:10]);% clabel(c,h);  
%   set(h,'EdgeColor','k') ; hold on
  
  hold on
 view(-0,90);
   
%    set(gca,'ylim',[ylm4(1)+0.015 ylm4(2)-0.015]);   set(gca,'xlim',[xlm4(1)+0.015 xlm4(2)-0.015]);
   set(gca,'ylim',[ylm4(1)-0.015 ylm4(2)+0.015]);   
   set(gca,'xlim',[xlm4(1)-0.015 xlm4(2)+0.015]);   
  xlim= get(gca,'xlim'); dx=(xlim(2)-xlim(1))*0.1; xt=xlim(1)+0.7*dx;
  ylim= get(gca,'ylim'); dy=(ylim(2)-ylim(1))*0.1; yt=ylim(1)-0.5*dy;
   
   ylim1=get(gca,'ylim'); dya=abs(diff(ylim1)*0.1)*0.25;
   xlim1=IM_range_lon;    ylim1=[IM_range_lat(1)-dya IM_range_lat(2)+0.01]  ;
set(gca,'xlim',xlim1,'ylim',ylim1);
 xlim= get(gca,'xlim'); dx=(xlim(2)-xlim(1))*0.1; xt=xlim(1)+0.7*dx;
 ylim= get(gca,'ylim'); dy=(ylim(2)-ylim(1))*0.1; yt=ylim(1)-0.5*dy;
                              cfm='%5.0f'; % cfm='%5.2f';
 ctexf =[ 'Sv ' num2str(mef ,cfm), ' \pm ' ,num2str(mvsf,cfm), ...
       ' [' num2str(mnf,cfm ) ',' num2str(mxf,cfm) '] ' ' '  ];  
%  text(xt,yt,ctexf,'FontSize',8);
hold on;
   
   set(gca,'TickDir','out'); %in
   set(gca,'YMinorTick','on','XMinorTick','on');
  daspect4 =[1  cosd(mean(ylm4)) 1];
  set(gca,'dataaspectratio',daspect4);
%     C2L =clabel(c,h,[6:1:16],'LabelSpacing',572,'FontSize',13);

  box on;
  caxis(caxD);
% if jp6==1,...     
    cbar=colorbar('location','EastOutside');  %cbposi = [0.83 0.13 0.015 0.640];  
                      cbrposi =[0.93086 0.24 0.015 0.570];
   if cm(1:3)=='zoo' ,cbrposi =[0.930 0.24 0.015 0.570]; end;
  set(cbar,'position',  cbrposi);
% end

  if do_log=='Y',...        
  %    lgwtck=[10.^(log10(10e-9):log10(-1))];
       lgwtck=[10.^(log10(ExpLim):log10(-1))];
    if jp6==1,...       
      set(cbar,'ytick',log10(lgwtck),'yticklabel',lgwtck,'tickdir','out');
                     cbt2=['C,log_{10}'];
   set(get(cbar,'title'),'string',cbt2);
    end
   caxis([log10(min(lgwtck)) log10(max(lgwtck))]);   
 end
   
%  colormap(jet);
 % colormap(vivid(12,'sr')); colormap(flipud(lbmap(11,'RedBlue')));
   colormap(flipud(lbmap(30,'Step5Seq')));   %salt75
%  colormap( (lbmap(30,'Step5Seq')));   
%  colormap(vivid); % ultraviolet oceandep for depth
   
   for i=1:size(l)-1
%    hfl= fill(coast(l(i)+1:l(i+1)-1,1),coast(l(i)+1:l(i+1)-1,2),[1 1 0.995], 'FaceAlpha', 0.25);
	      fill(coast(l(i)+1:l(i+1)-1,1),coast(l(i)+1:l(i+1)-1,2),[1 1 1], 'FaceAlpha', 0.5);
   end
   
  %% google-maps use here
    box on;
       cbfvw=datestr((FVCOM1.Date(it)),' yyyymmdd-HH');
    dx=(xlm4(2)-xlm4(1))/10; dy=(ylm4(2)-ylm4(1))/10;

  clear uvm mnu mxu meu mvr uvm1 Cax zz
%                  zz=varn;
                 zz=varmp;
 uvm=zz(isfinite(zz));
 mnu=min( nanmin(uvm)' ); mxu=max( nanmax(uvm)' ); 
 meu=mean(nanmean(uvm)'); mvr=var(uvm); mvs=std(uvm);  
 
   c_units =['']; cvr=['Sv']; cfm='%5.0f'; % cfm='%5.3g';
   ctexv =[  num2str(meu,cfm), ' \pm ' ,num2str(mvs,cfm), ...
        ' [' num2str(mnu,cfm) ',' num2str(mxu,cfm) '] ' c_units  ];  
%   plot(MOR.lon,MOR.lat,'+k','markersize',1.5); hold on
  
  try ,
      plot(GA.lon,GA.lat,'ok','markersize',1.2); hold on
%     plot(CTD.lon,CTD.lat,'ok','markersize',1.2); hold on
%     text(CTD.lon,CTD.lat,CTD.ID); hold on
  end
 
   title_txt=['Sv at iz=' ciz ',' cbfvw  ' ' ];
%  set(htl,'position',[ xlm4(1)+3.25*dx   ylm4(2)+0.01  1.0])  
%      text(xlm(1)+7.0*dx, ylm(2)+0.99*dy,ctexv,'FontSize',9,'BackgroundColor',[1 1 0.98]); hold on
   htc=text(xlm4(1)+6.05*dx ,ylm4(2)-0.16*dy,ctexf,'FontSize',6 ,'BackgroundColor',[1 1 0.98]); hold on
   htl=text(xlm4(1)+0.01*dx ,ylm4(2)-0.21*dy,title_txt,'FontSize',7,'BackgroundColor',[1 1 0.98]); hold on

%% ------------------------   incut winds ------
if 0,...
 %   w_elm=Pr.cel(35); %Wind cell location;
 %   w_elm=3304; w_nod=1840;
  ax1=gca;
  posw=get(gca,'Position'); 
          posw3=[pos6(1)+0.165 pos6(2)+0.075 0.10 0.12];
%  if cm(1:3)=='zoo' ,    posw3=[0.73 0.32 0.10 0.12]; end
  ax3 = axes('Position',posw3);
clear uw vw uvw;
          uw=FVCOM1.uwind_speed(it,w_elm); 
          vw=FVCOM1.vwind_speed(it,w_elm);
          uvw=sqrt(uw*uw + vw*vw);
hpw=quiver(0,0, uw, vw,'sc','filled');
set(ax3,'ylim',[-4 4]);set(ax3,'xlim',[-4 4]);
set(hpw,'color','m','linewidth',2,'Markersize',3,'MarkerEdgeColor','c');
htw=text(-1.415,-2.095, ['' num2str(uvw,'%3.1f ') ' ms^{-1}'],...    
    'FontSize',11,'color','k','FontWeight','bold'); %ba

set(ax3,'visible','off'); %,'on'); 
set(ax3,'Position',posw3);
end
%% ==============     
 ax1=gca; 
%    dya=0.005; %; dya=abs(diff(ylim1)*0.1)*0.25;
     dya=0.01; %; dya=abs(diff(ylim1)*0.1)*0.25;
     ylm=[floor(ylim(1)-dya) ,ceil(ylim(2)+dya)];
     xlm=[floor(xlim(1)-dya) ,ceil(xlim(2)+dya)];
  set(ax1,'xtick',[xlm(1):0.1:xlm(end)]);
  set(ax1,'ytick',[ylm(1):0.1:ylm(end)]);
  ax1.XRuler.MinorTickValues=[xlm(1):0.05:xlm(end)];
  ax1.YRuler.MinorTickValues=[ylm(1):0.05:ylm(end)];
  ax1.XRuler.MinorTickValuesMode ='manual';
  ax1.YRuler.MinorTickValuesMode ='manual';
  
 set(gcf,'color','w'); % 'none')
 orient portrait ;
 set(gcf,'PaperPositionMode','auto');
 drawnow    
 fnm=([path_fig 'FV-Sv-iz' ciz '_' datestr(FVCOM1.Date(it),'YYYYmmdd_HH') '_sh' ]);    
 if do_log=='Y';  fnm= [fnm 'L' ];end;
%% %%%%%%%%%%%making gif vect:
  if  mkgif=='Y', ....
   if gi==0,...
%  nameuv=strcat([ FVCOM_plot_dirg 'vec_' fig_name var_plot cm ] );
          date_daye=[datestr(Date(1),'yymmdd') '_' datestr(Date(end),'mmdd')];
   nameu=strcat([ path_fig  'Sv_iz-' ciz  '_' date_daye    cm ] );
   if do_log=='Y';  nameu = [nameu  'L' ];end;   
                   gifname=[nameu '.gif'];
   II2 = getframe(gcf);
   II2 = frame2im(II2);
   [XX2, map] = rgb2ind(II2, 256); % need  Image Processing Toolbox
   imwrite(XX2, map, gifname, 'GIF', 'WriteMode', 'overwrite', 'DelayTime', 0, 'LoopCount', Inf);
   gi=gi+1;
   else
   II2 = getframe(gcf);
   II2 = frame2im(II2);
   [XX2, map] = rgb2ind(II2, 128);%
   imwrite(XX2, map, gifname, 'GIF', 'WriteMode', 'append', 'DelayTime', 0);
   gi=gi+1;
   end; %gi
   clear II2 XX2
  end; %mkgif
%% =================  
    fnm=([path_fig 'Sv-iz' ciz '_' datestr(FVCOM1.Date(it),'YYYYmmdd_HH') '_sh' cm ]);     
 if jp6==6, jp6=0;    
   if ispc,
   set(gcf,'renderer','zbuffer') ;  
%  set(gcf,'renderer','opengl')
   else
      set(gcf,'renderer','painters') ;        
   end  
      print(['-f'],'-dpng','-r300',[fnm  '.png']);     %   700k
%       stop
 end; %jp6    
   if jp6<=2,   print(['-f'],'-dpng','-r300',[fnm  '.png']);  end %   700k

end

disp([fnm]);
end
