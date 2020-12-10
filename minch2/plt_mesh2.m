 figure(2);clf
 addpath(  'M:/Mar_Phys/matlab/m_map/' );
 FVCOM_in_dir =['../Archive/mat_2016/'];

 load([ 'Mesh.mat'])
  [IM_range_lat ,IM_range_lon] = deal([54.45  58.70],[ -7.65   -4.55]);
 m_proj('Lambert Conformal Conic','lon',[ IM_range_lon],'lat',[IM_range_lat],      'ell','wgs84')
 m_grid('box','fancy')
%plot vertices
   col=0.7*[1 1 1];
   [X,Y]=m_ll2xy(Mesh.geog(:,1),Mesh.geog(:,2),'clip','on');

   jnk=patch('Vertices',[X,Y],'Faces',Mesh.trinodes,...
	   'EdgeColor',col,'FaceColor','none');
hold on
% m_plot(Probe.nodexy(:,1), Probe.nodexy(:,2),'sm'...
%     ,'LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','g','MarkerSize',7);
% hold on;
% patch('Vertices',[Mesh.nodexy(:,1),Mesh.nodexy(:,2)],'Faces',Mesh.trinodes,...
%        'Cdata',Mesh.depth(:,3),'edgecolor','interp','facecolor','interp');

set(gcf,'position',[53  53  480   850]);
view(0,90);
% m_text(Probe.nodexy(:,1)+0.05, Probe.nodexy(:,2)+0.05, char(Probe.name(:)) )
set(gcf,'PaperPositionMode','auto')
fig_name=['minch2_Mesh'];
print(['-f2'],'-dpng','-loose','-r500',...
[ 'Map_' fig_name   '.png']);
