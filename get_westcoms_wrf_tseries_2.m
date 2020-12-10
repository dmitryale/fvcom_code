% get_westcoms_wrf_tseries_2
% 1. download the several latest WeStCOMS-WRF file(s) from SAMS Thredds server:
% 2. identify target location x0,y0 ( "meteo station" )
% 3. derive data for those site(s) from nearest WRF grid points
% PS 4.  do not use experimental fields : cvars_rm={'sst';'tcc';'cld_fra';'qvapor'};
%# ./wget.exe http://sausage.sams.local:8081/thredds/fileServer/scoats/WRF/Archive/netcdf_2019/wrf_20191023-1030_d03.nc
%             http://sausage.sams.local:8081/thredds/catalog/scoats/WRF/Archive/netcdf_2019/catalog.html
%(c)dmitry.aleynik@sams.ac.uk: 5 November 2019
%% ><(((((({°> ~~~~~.~~~.~
% Copyright (c) 2019, SAMS ;#C:\samhanach\sa01da\work\WRF\WRFout\OUT
%======================================================================
  close all; clear all;
  this_dir=pwd; DRV=this_dir(1:1);
  addpath(['M:/Mar_Phys/matlab/general/']);

   cvars_rm={'sst';'tcc';'cld_fra';'qvapor'}; % should not be used NB
   those_dir=['W:/sa01da/work/WRF/Archive/'];
   dir_wrf_nc=dir([ this_dir '/netcdf-2019/w*_d03.nc']);% ->20190612
%   dir_wrf_nc=dir([ this_dir '/netcdf_2019/w*_d03.nc']);
[nf,fi]=size(dir_wrf_nc);
 for fi=1:nf
     fname=[dir_wrf_nc(fi).folder '/' dir_wrf_nc(fi).name];
        [finfo,outstrct ] = read_nc_file_struct(fname);
try               outstrct=rmfield(outstrct,cvars_rm); end % removing 'experimental' fields
  A(fi,1).WRF=outstrct;
  end
    time_offset = 678942; % from FVCOM time to matlab time

 %% define station location to derive meteo data:
%st.x(1)=-5-26.417/60;  st.y(1)=56+26.986/60; st.name(1,:)=['Duns']; % DUSNTAFFNAGE meteo
%  st.x(1)=-5-30.12/60;   st.y(1)=56+28.73/60;  st.name(1,:)=['LY-1']; % Lismore sta LY1
%st.x(2)=-5-29.00/60;   st.y(2)=56+25.00/60;  st.name(2,:)=['Oban']; % Oban tides 5+29/60
%st.x(3)=-6-26.00/60;   st.y(3)=56+35.00/60;  st.name(3,:)=['Tire']; % Tiree mooring, dep 45 m
%st.x(4)=-5-26.1865/60; st.y(4)=56+27.208/60; st.name(4,:)=['Pont']; % DUSNTAFFNAGE Pontoon

%st.x(1)=-5.434300;   st.y(1)=56.480567 ;  st.name(1,:)=['AB-1']; st.name_long(1,:)  ='Ardmucknish Bay';
 st.x(1)=-5.329133;   st.y(1)=56.622217 ;  st.name(1,:)=['BG-2']; st.name_long(1,:)='Balnagowan     ';
%	Ardmucknish Bay -5.4343   56.480567   ; # 3/4/2019 to 19/4/2019     &     6/5/2019 to 24/5/2019
% Balnagowan      -5.329133 56.622217   ; # 14/10/2019 to 29/10/2019
%#-6.272500, 57.045000, 31.0, 20190329-00:00, 20200817-00:00, NE-Rum         , 20190513 01:14, 07:14, t,s,u,v,el,uw,vw,29274, 51673
%#-5.586017, 56.185866, 20.0, 20190329-10:53, 20191114-00:00, BDNC           , 20190513 01:14, 07:14, t,s,u,v,el,uw,vw, 1419,  2548
       x0=st.x'; y0=st.y';
 WRFst.x0=x0;
 WRFst.y0=y0; WRFst.sta = st.name;
   [index,distance]=nearxy(A(1).WRF.XLONG(:),A(1).WRF.XLAT(:),x0,y0); %index=483
            [i0,j0]=ind2sub(size(A(1).WRF.XLONG),index);
 WRFst.index=[i0,j0,index];
%% verify the indexes for given site, scater-interpolant or splines would be better then nearest grid point,  but a bit slower
  if 0
    figure(1);clf;
    plot(A(1).WRF.XLONG,A(1).WRF.XLAT,'.'); hold on
    plot(A(1).WRF.XLONG(i0,j0),A(1).WRF.XLAT(i0,j0),'pm'); hold on
    plot( x0, y0,'rs'); hold on
    load('M:\Mar_Phys\matlab\m_map\uk_coastline.mat');  l=find(isnan(coast(:,1)));
  for i=1:size(l)-1
%    hfl= fill(coast(l(i)+1:l(i+1)-1,1),coast(l(i)+1:l(i+1)-1,2),[1 1 0.995], 'FaceAlpha', 0.25);
	      fill(coast(l(i)+1:l(i+1)-1,1),coast(l(i)+1:l(i+1)-1,2),[1 1 1], 'FaceAlpha', 0.5); hold on
  end
          namepng0 =[ 'map_WRF_' st.name   ];
print(['-f'],'-dpng','-loose','-r500',[namepng0  ,'.png']) ;

  end
%% ---------------------

 WRFst.Times(1,:)=A(1).WRF.Times(:,1)';
 WRFst.Times(2,:)=A(1).WRF.Times(:,2)';
 WRFst.Times(3,:)=A(1).WRF.Times(:,end)';

     cvars   = fieldnames(A(1).WRF); % cvars=cvars(5:end);
 for fi=1:nf
       mtime=datenum(A(fi).WRF.Times');
  go=[]; a=[];
    if fi==1
       WRFst.mtime=mtime; go=[1:length(mtime)]';
    else
        go=find(mtime > WRFst.mtime(end));
    end
       if isempty(go), continue;end
    if fi >1,  WRFst.mtime=cat(1,WRFst.mtime, mtime(go)); end
   for iv=1:length(cvars)
       cvar=char(cvars(iv));
                 a=double(A(fi).WRF.(cvar));
       dims=size(a);
       if length(dims)<3, continue ;end
       if fi==1
    WRFst.(cvar)=(squeeze(a(i0,j0,go)));
       else
    WRFst.(cvar)=cat(1,WRFst.(cvar),squeeze(a(i0,j0,go)));
       end
   end
 end

    db=datestr(WRFst.mtime(1),'yyyymmdd');
    de=datestr(WRFst.mtime(end),'mmdd');
save(['tmseries_WRF_' st.name '_' db '-' de '.mat' ],'WRFst','-v7.3');

 if 1
 figure(2);clf; set(gcf,'position',[50 200 800 250]);
                                scaleX=cosd(mean(y0));
                                yy=zeros(length(WRFst.U10),1);
 quiver( WRFst.mtime,yy,  WRFst.U10,  WRFst.V10*scaleX, 1 );  hold on;
 axis equal; datetick('x'); axis tight;  box off; set(gca,'Ytick',[ ])
  set(gcf,'color','w'); % set(gcf,'color','none')
title (['WRF UV10m at ' st.name ]);
          namepng =[ 'tmseries_WRF_' st.name '_UV10_' db '-' de ];
print(['-f'],'-dpng','-loose','-r500',[namepng  ,'.png']) ;

 end
