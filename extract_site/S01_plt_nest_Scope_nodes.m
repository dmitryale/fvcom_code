% function[]= S01_plt_nest_Scope_nodes()
% INput :  input_lag_Scope.csv 
% OUTput:  minch2_Scope_nest.mat, minch2_Scope_nestnodes.dat etc 
% based on plt_nest_Scope
% created : Dmitry Aleynik , 29/10/2018  
%  __o_O__/X/                                                             &
%  \_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS, COMPASS/CAMPUS 2018-2021     &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%=========NCNEST
clear all; 
  a=pwd;  DRV=a(1:1); this_dir = a; 
  addpath(['M:/mar_phys/matlab/seawater/']);
  addpath('../');
  
             data_dir =[ this_dir '/' ] ; % '/Scope/'];
 filename =[ data_dir '/input_lag_Scope.csv'];   
 clearvars a OU OUT Pr Pr0 
 
 mk_fig='Y';
 OU = csvread(filename);
 OUT.lon=OU(:,1);
OUT.lat=OU(:,2);
OUT.dep=OU(:,3);
  xlmd=[-7.0 -4.90]; ylmd=[55.7 56.70];
% load('..\Scope\Jura_20160728.mat') ; mtime=datenum(JuraRawData.DateandTimeGMT);
                    cnm1='Loch-Buie';
                    cnm2='MenardPoint';
                    cnm3='Jura-West';
for ii=1:length(OU(:,1)),...
        a=[]; a=sprintf('%d',OU(ii,4));
        b=[]; b=sprintf('%d',OU(ii,5));
        if ii==1, cnm=cnm1;end
        if ii==2, cnm=cnm2;end
        if ii==3, cnm=cnm3;end
OUT.ID(ii,:) =[ cnm(1:6) a];
OUT.mtb(ii,:) =datenum(a,'yyyymmdd');
OUT.mte(ii,:) =datenum(b,'yyyymmdd');
end
OUT.Site=OUT.ID;

 bn(1,:)=[ 1,5,10];
 Pr.lon=OUT.lon';   Pr.lat=OUT.lat';
 Pr.x  =OUT.lon';   Pr.y  =OUT.lat';
 Pr.dep=OUT.dep;    Pr.ID =OUT.ID;
 Pr.name=cnm;       Pr.bn = bn;
 Pr.mtb=OUT.mtb;    Pr.mte=OUT.mte;
% load minch2 mesh:
  Minch2=load([  'c:/samhanach/sa01da/work/minch2/code/mesh.mat']);
  try Minch2.mesh.tri=Minch2.mesh.trinodes;end
  mesh=Minch2.mesh;
   tri=mesh.tri;
     mesh.x=mesh.geog;
     mesh.x(:,3)=mesh.depth;

%% ===================
 data_inp=data_dir; % 'in/';
 path_fig=data_inp;

if 1,...
  for i=1:length(Pr.x),...
          x0(i,1)=Pr.x(i);
          y0(i,1)=Pr.y(i);
    clear rr jc rmin dd xx yy zz ;
   for j=1:mesh.Nverts,...
           clear dd xx yy ;
        xx=mesh.x(j,1);
        yy=mesh.x(j,2);
%       zz=mesh.x(j,3);
                    dd=sqrt( ((xx-x0(i))^2) +  ((yy-y0(i))^2) );
            rr(j,1)=dd;
                                 la=[mesh.geog(j,2), Pr.lat(i)];
                                 lo=[mesh.geog(j,1), Pr.lon(i)];
      [dist,phaseangle]=sw_dist(la,lo,'km');
            rd(j,1)=dist;
   end

       [rmin, jc] =  min(rr);
     % [rMin jC]  =  sort(rr);
       [rMin, jC] =  sort(rd);
% Get 3 nearest
           nj=[1:6]';
%          nj=[1:5]';
%          nj=[1:3]'; % nearest 3

%   Pr.jc(i,1)=jc(1);
    Pr.jc(i,:)=jC(nj);
    Pr.xx(i,:)=mesh.x(jC(nj),1);
    Pr.yy(i,:)=mesh.x(jC(nj),2);
    Pr.zz(i,:)=mesh.x(jC(nj),3);
    Pr.rd(i,:)=rMin(nj);
    Pr.rr(i,:)=rmin(1 );
    Pr.id_3(i,:)=repmat(i,[length(nj),1])';
   % Pr.difz(i,:)=min( abs(Pr.zz(i,:)-Pr.dep(i)));

  end

%==========
% get extract unique only
Pr0=Pr;

[ujc,ujci]=unique(Pr0.jc);
%  exclude one maximum
   rd=Pr0.rd(ujci); [rdmx, excl]=max(rd);
      a=ujci(excl);
      ujci=ujci(ujci~=a);

Pr.jc=Pr0.jc(ujci);
Pr.xx=Pr0.xx(ujci); Pr.yy=Pr0.yy(ujci);
Pr.zz=Pr0.zz(ujci); Pr.rd=Pr0.rd(ujci);
Pr.id=Pr0.id_3(ujci);
% Pr.rr=Pr0.rr(ujci);


x=mesh.x(:,1);
y=mesh.x(:,2);
          dims = size(mesh.tri);

   mark = ones(length(x),1);

clear rd jC jc
  for i=1:length(Pr.x),...
          x0(i,1)=Pr.x(i);
          y0(i,1)=Pr.y(i);
    clear rr jc rmin dd xx yy zz ;
   for j=1:mesh.Nelems,...
           clear dd xx yy lo la ;
        xx=mesh.uvnode(j,1);
        yy=mesh.uvnode(j,2);
%       zz=mesh.uvdepth(j,3);
                    dd=sqrt( ((xx-x0(i))^2) +  ((yy-y0(i))^2) );
            rr(j,1)=dd;
                                 la=[mesh.uvnode(j,2), Pr.lat(i)];
                                 lo=[mesh.uvnode(j,1), Pr.lon(i)];
      [dist,phaseangle]=sw_dist(la,lo,'km');
            rd(j,1)=dist;
   end

        [rmin, jc] =  min(rr); %degrees
        [rMin, jC] =  sort(rd);
% Get 3 nearest
%              nj=[1:3]';
                nj=[1:6]'; % nearest 6
%               nj=[1:5]'; % nearest 5
 %              nj=[1:3]'; % nearest 3
     Pr.jce(i,:)=jC(nj);
     Pr.rdc(i,:)=rMin(nj);
     Pr.idc_3(i,:)=repmat(i,[length(nj),1])';
   end

  Pr0.cel=Pr.jce;
  Pr0.jce=Pr.jce;
  Pr0.idc_3=Pr.idc_3;

[uce,ucei]=unique(Pr0.cel);
Pr.cel =Pr0.cel(ucei);
Pr.jce =Pr0.jce(ucei);
Pr.idce=Pr0.idc_3(ucei);

%chek that  noj=Pr.nod:
%   noj=find(mark==0);

    Pr.nod=Pr.jc;
  % dump a probe node & cell file to ascii
  numN=prod(size(Pr.nod));
  numC=prod(size(Pr.cel));
num=max(numC,numN);

  for ii=1:length(Pr.xx)
      if ii<=length(Pr.cel),
  plot(mesh.uvnode(Pr.cel(ii),1),mesh.uvnode(Pr.cel(ii),2), ...
      'r^','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',8.5);hold on
      end
  plot(Pr.xx(ii),Pr.yy(ii), ...
      'b.','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6.5);hold on
  if ii<=length(Pr.lon),
  plot(Pr.lon(ii),Pr.lat(ii), ...
      'k+','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8.5);hold on
  text(Pr.lon(ii)+0.001,Pr.lat(ii)+0.001, ...
      [ num2str(ii,'%3.2d') ] );hold on
  end
  end


if 1,
     fid = fopen([data_inp 'probes' num2str(length(Pr.nod)) '_nml_stnames.csv'],'w');
     fprintf(fid,'%s, \n', [ 'N#_mesh, Node_id, Cell_id, N#_site, Site_name'] );
for i=1:num,
    if i<=numC,cei=Pr.cel(i);end
    if i<=numN,cni=Pr.nod(i);end
       PR_ID=[]; PR_ID=Pr.idce(i);
       site_name =[]; % site_name = char( OUT.Site{PR_ID,1} ) ;
       site_name = OUT.Site(PR_ID,:) ;
  fprintf(fid,'%3d, %6d, %6d, %4d, %s, \n', i ,cni, cei, PR_ID, site_name );
end;
fclose(fid);
% stop;
end

%%%
     fid = fopen([data_inp 'probes' num2str(length(Pr.nod)) '_nml.dat'],'w');
for i=1:num,
    if i<=numC,cei=Pr.cel(i);end
    if i<=numN,cni=Pr.nod(i);end
  fprintf(fid,'%3d %6d %6d \n', i ,cni, cei );
end;
fclose(fid);

% dump_nest files
     fid = fopen([data_inp 'minch2_Scope_nestnodes.dat'],'w');
  fprintf(fid,'%s %3d \n','Node_Nest Number =', numN );
 for i=1:numN
  fprintf(fid,'%3d %6d %3d \n', i ,Pr.nod(i), i );
end;
fclose(fid);

     fid = fopen([data_inp 'minch2_Scope_nestcells.dat'],'w');
  fprintf(fid,'%s %3d \n','Cell_Nest Number =', numC );
 for i=1:numC
  fprintf(fid,'%3d %6d %3d \n', i ,Pr.cel(i), i );
end;
fclose(fid);


dir_nestc=dir([data_inp '/mi*nestcells.dat']);
           fnc=deblank(dir_nestc(1,:).name);
 newData1 = importdata([data_dir fnc], ' ',1 );
 nest.CID_mi(:,1)=newData1.data(:,2); % minch2 cells
 nest.CID(:   ,1)= nest.CID_mi(:,1); %  obc_elems;

  dir_nest=dir([data_inp '/mi*nestnodes.dat']);
            fnn=deblank(dir_nest(1,:).name);
  newData1 = importdata([data_dir fnn], ' ',1 );
 nest.NID_mi(:,1)=newData1.data(:,2); % minch2 cells
 nest.NID(:   ,1)= nest.NID_mi(:,1); %  obc_elems;
[nl, ni ]= size(newData1.data);
nl=[1:nl];
  nest.mtb=Pr.mtb;
  nest.mte=Pr.mte;

  [obc_nodes] = nest.NID_mi;
  [obc_elems] = nest.CID_mi;

% [nest.num, a]=size(dir_nest);
 nest.num=1;
for mi=1:nest.num,...
%          clear data textdata fileToRead1 newData1 nl  vars
nest.np( mi,1 )   = max(nl);
nest.x(  nl,:)    = mesh.x(   nest.NID(nl,mi),1:3);
nest.xg( nl,:)    = mesh.geog(nest.NID(nl,mi),1:2);

              ce=tri(obc_elems,:);
   ncl=length(ce(:,1));
   for ic=1:ncl,
   nest.xgc(ic,1)   = mean(mesh.geog(ce(ic,:),1));
   nest.xgc(ic,2)   = mean(mesh.geog(ce(ic,:),2));
   end

% external nodes, cels:
   nest.xg_mi(nl,:) = Minch2.mesh.geog(nest.NID_mi(nl,mi),1:2);
   cex=Minch2.mesh.tri(nest.CID_mi,:); ncli=length(cex(:,1));
   for ic=1:ncli,
   nest.xgc_mi(ic,1)   = mean(Minch2.mesh.geog(cex(ic,:),1));
   nest.xgc_mi(ic,2)   = mean(Minch2.mesh.geog(cex(ic,:),2));

%    text(nest.xgc_mi(ic,1),nest.xgc_mi(ic,2), ...
%       {(nest.CID_mi(ic,1))} );hold on

   end
end


 clear data textdata fileToRead1 newData1 nl vars namepng
   nest.obc_nodes=obc_nodes;
   nest.obc_elems=obc_elems;
%% ==========
if mk_fig=='Y',...
  col07=0.7*[1 1 1];
  col09=0.9*[1 1 1];
%
figure (125); clf
% set(gcf,'Position',[55 55 700 700]);
  set(gcf,'Position',[55 55 450 600]);

% patch('Vertices',[mesh.geog(:,1),mesh.geog(:,2)],'Faces',mesh.tri,...
%           'Cdata',-mesh.x(:,3),'edgecolor','interp','facecolor','interp'); hold on


%plot(nest.x(:,1),nest.x(:,2), ...
 plot(nest.xg(:,1),nest.xg(:,2), ...
      'kd','MarkerEdgeColor','r','MarkerFaceColor','m','MarkerSize',4.5);hold on
 plot(nest.xg_mi(:,1),nest.xg_mi(:,2), ...
      'gs','MarkerEdgeColor','y','MarkerFaceColor','c','MarkerSize',3.5);hold on

 plot(nest.xgc(:,1),nest.xgc(:,2), ...
      'kv','MarkerEdgeColor','y','MarkerFaceColor','c','MarkerSize',3.5);hold on
 plot(nest.xgc_mi(:,1),nest.xgc_mi(:,2), ...
      'k^','MarkerEdgeColor','c','MarkerFaceColor','y','MarkerSize',4.5);hold on

  for ii=1:length(Pr.x)
  plot(Pr.lon(ii),Pr.lat(ii), ...
      'k+','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8.5);hold on
  text(Pr.lon(ii)+0.01,Pr.lat(ii)+0.01, ...
      [ num2str(ii,'%3.2d') ] );hold on
  end

axis equal
% colorbar
title('NCNEST_ON  large 2 small domain : nodes/elems ')

 xlmg=[min(mesh.geog(:,1)) max(mesh.geog(:,1)) ];
 ylmg=[min(mesh.geog(:,2)) max(mesh.geog(:,2)) ];

 xlm=[min(mesh.x(:,1)) max(mesh.x(:,1)) ];
 ylm=[min(mesh.x(:,2)) max(mesh.x(:,2)) ];
%  set(gca,'xlim',xlm); set(gca,'ylim',ylm)
 set(gca,'xlim',xlmg); set(gca,'ylim',ylmg)
 daspect([1 cosd(mean(ylmg)) 1]);
% set(gcf,'position',[50   100   525  900]);
 cbar=colorbar;
 cbp=[0.91 0.34    0.03    0.25]; set(cbar,'position',cbp);

box on;
  set(gca,'YMinorTick','on','XMinorTick','on');
  set(gca,'TickDir','out'); %in
  set(gcf,'PaperPositionMode','auto');
jnk=patch('Vertices',[mesh.geog(:,1),mesh.geog(:,2)],'Faces',mesh.tri,...
	      'EdgeColor',col07,'FaceColor','none'); hold on;
 namepng=[path_fig 'NCnest_fileg']; % saveas(gcf,namepng,'png') ;
 print(['-f'],'-dpng','-loose','-r300',[namepng  ,'.png']) ;
 
 set(gca,'xlim',xlmd);   set(gca,'ylim',ylmd);
 namepng=[path_fig 'NCnest_filegZ']; % saveas(gcf,namepng,'png') ;
 print(['-f'],'-dpng','-loose','-r300',[namepng  ,'.png']) ;

%% -----zoom it

for ii=1:length(Pr.x),...

figure (126); clf
    set(gcf,'Position',[55 55 600 450]);
 plot(nest.xg(:,1),nest.xg(:,2), ...
      'kd','MarkerEdgeColor','r','MarkerFaceColor','m','MarkerSize',4.5);hold on
 plot(nest.xg_mi(:,1),nest.xg_mi(:,2), ...
      'gs','MarkerEdgeColor','y','MarkerFaceColor','c','MarkerSize',3.5);hold on

 plot(nest.xgc(:,1),nest.xgc(:,2), ...
      'kv','MarkerEdgeColor','y','MarkerFaceColor','c','MarkerSize',3.5);hold on
 plot(nest.xgc_mi(:,1),nest.xgc_mi(:,2), ...
      'k^','MarkerEdgeColor','c','MarkerFaceColor','y','MarkerSize',4.5);hold on

xlmgz=[min(Pr.x(ii))-0.0525  max(Pr.x(ii))+0.0525];
ylmgz=[min(Pr.y(ii))-0.025  max(Pr.y(ii))+0.025 ];

% xlmgz=[min(nest.xg(:,1))-0.125 max(nest.xg(:,1))+0.125];
% ylmgz=[min(nest.xg(:,2))-0.05 max(nest.xg(:,2))+0.05];
set(gca,'xlim',xlmgz);   set(gca,'ylim',ylmgz);
jnk=patch('Vertices',[mesh.geog(:,1),mesh.geog(:,2)],'Faces',mesh.tri,...
	      'EdgeColor',col07,'FaceColor','none'); hold on;

%    mk_topo='N';
   mk_topo='Y';
 if ii==1 & mk_topo=='Y',...
      vrmf= -mesh.depth;
% vs10=[-250:10:0  ];
  vs05=[-250:05:0  ];
%   mnfz =min( nanmin(vrmf)' );    mxfz  =nanmax( nanmax(vrmf)' );
%   menz =nanmean(nanmean(vrmf)'); mvsz  =nanmean(std(vrmf));
  end;
 [c2,h2] = tricontour(mesh.geog,mesh.trinodes,vrmf,vs05); hold on;
  set(h2,'LineWidth',0.10,'linestyle','-');%,'edgecolor','k');col09); %...'k');
%[c3,h3] = contour(mesh.geog(:,1),mesh.geog(:,2),vrmf,vs10);hold on; % shading interp; hold on;
  caxis([-250 0]);
 cbar=colorbar;
 cbp=[0.88 0.20 0.02 0.35]; set(cbar,'position',cbp);
 cbar.Label.String=['Depth,m'];

  plot(Pr.lon(ii),Pr.lat(ii), ...
      'k+','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8.5);hold on
  text(Pr.lon(ii)+0.001,Pr.lat(ii)+0.001, ...
      [ num2str(ii,'%3.2d') ] );hold on

     daspect([1 cosd(mean(ylmg)) 1]);

set(gcf,'PaperPositionMode','auto');
 namepngi=[path_fig 'NCnest_fileg_' num2str(ii,'%2.2d') ];
print(['-f'],'-dpng','-r300',[namepngi '_zl' '.png']) ;
end

end
  save(['minch2_Scope_nest.mat'],  'nest');
else
  load(['minch2_Scope_nest.mat']);% ,  'nest');
end
  Pr.nod = nest.obc_nodes;
  Pr.cel = nest.obc_elems;
%% ==================
