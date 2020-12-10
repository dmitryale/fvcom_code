function [FVCOM1]=get_FVCOM1_fnc(fnam )
%Usage example
%     addpath (['../extLib/']);
%    FVCOM_inc_dir =['../Archive/netcdf_2016/'];
%    dif_Fin =  dir([FVCOM_inc_dir 'm*2_2*.nc']);
%    nf=length(dif_Fin); f1=1;
%           f1=nf; % get the last one
%    for fi=f1:nf
%       fnam = [FVCOM_inc_dir dif_Fin(fi,1).name];
% %      [FVCOM1]=get_FVCOM1_fnc(fnam );
% [finfo outstrct]=read_nc_file_struct(fnam);
%    end
%   load([ 'Mesh.mat'])
%(c)dmitry.aleynik@sams.ac.uk, 2015.09.14)
%% ><(((((({°> ~~~~~.~~~.~
global Mesh
  time_offset = 678942; % from FVCOM time to matlab time
  slayer='tanh'; caseN ='2'; casename=['minch' caseN];
  in_cdf=1;

 if ispc,... % isunix
  addpath (['M:/Mar_Phys/matlab/seawater/']);
 else
  addpath('/home/sa01da/Mar_Phys/matlab/seawater/');
 end
 
  addpath (['../extLib/']);
 
[finfo outstrct]=read_nc_file_struct(fnam);
W = outstrct; clear outstrct;
    if W.lon(1)>180,...
     W.lon=W.lon-360;
     W.lonc=W.lonc-360;
    end
      W.xc=W.lon;%c;
      W.yc=W.lat;%c;
     date_s=W.Times'; 
   [nt sti]=size(date_s);


 FVCOM=W; FVCOM.el=FVCOM.zeta; FVCOM.t=W.temp;FVCOM.s=W.salinity;
 
 Date = double( FVCOM.Itime2 ./ (3600000) )./24 + double(FVCOM.Itime) + time_offset ;
 nt=length(Date);
 FVCOM.siglay=W.siglay; %10
 FVCOM.siglev=W.siglev; %11
     ZZ=FVCOM.siglay(:,:);
  [m nz]=size(FVCOM.siglay);
  
 clear  FVCOM_data  RW W ;
   [i1 i2 i3]=size(FVCOM.('t'));
                   ord2d=[2,1];
      if i1<i3,    ord3d=[1,2,3]; end; %keep as it is
      if nt==1 || i1>i3,
                   ord3d=[3,2,1];
      FVCOM.('el')=permute(FVCOM.('el') ,ord2d);
      end

 if  nt==1,% singleton
 el= FVCOM.el(:)';
 else
 el= FVCOM.el(:,:);
 end

for it=1:nt,...
%z=FVCOM.h(:)+FVCOM.zeta(:,tt);
 if  nt==1,% singleton
 z  =FVCOM.h(:,1)+el(:);
 else
 z  =FVCOM.h(:,1)+el(it,:)';
 end

 for iz=1:nz,...
     if isempty(slayer)~=1 && strcmp(slayer(1:4),'unif')==1,...
       if iz==1,...
       FVCOM.z(it,iz,:)= - el(it,:)';
       else
       FVCOM.z(it,iz,:)=iz.*(z./sigma_levs)-(z./sigma_levs)*0.5;
       end
     end;

     if isempty(slayer)~=1 && strcmp(slayer(1:4),'tanh')==1,...
       if iz==1,...
       FVCOM.z(it,iz,:)=  - el(it,:)';
       else
       FVCOM.z(it,iz,:) = -ZZ(:,iz).*z;         % depth on sigma level nodes
       end

       if (iz==nz)    ,...
       FVCOM.z(it,iz,:) = FVCOM.h(:,1);            % bottom depth
       end

     end;

end
end;
        ntime0 = length(FVCOM.time(:,1));
      FVCOM.('t') =permute(FVCOM.('t')  ,ord3d);
      FVCOM.('s') =permute(FVCOM.('s')  ,ord3d);

  for it=1: ntime0,...
  clear rr zz ss tt pres depth mlat;

  tt=[];tt(:,:)=squeeze(FVCOM.('t')(it,:,: )); tt=double(tt)';
  ss=[];ss(:,:)=squeeze(FVCOM.('s')(it,:,: )); ss=double(ss)';
  zz=[];zz(:,:)=squeeze(FVCOM.('z')(it,:,: )); zz=double(zz)' ;

%  ntr =  length(tt(:,1));rlat=mean(Mesh.geog(trn_nodes,2)); rlat=mean(Mesh.geog(:,2);
   rlat=Mesh.geog(:,2); ntr=1;
  mlat=[]; mlat(:,:)=repmat(rlat,ntr,nz);
  depth=zz;
             pres(:,:)= sw_pres(depth(:,:),mlat) ;
  rr=[];rr(:,:)=sw_pden(ss,tt,pres,0)-1000.0; %pot_sigma
  FVCOM.r(it,:,:)=single(rr(:,:) )';
  end; %it  %tanh

   if in_cdf==1  || out_mat==1 ,...  ,
    FVCOM1.filename = FVCOM.filename; %flnm
    FVCOM1.Date   = Date;
    FVCOM1.time   = FVCOM.time  ;
    FVCOM1.Itime  = FVCOM.Itime ;
    FVCOM1.Itime2 = FVCOM.Itime2 ;
    FVCOM1.s = permute(FVCOM.s, ord3d );
    FVCOM1.t = permute(FVCOM.t, ord3d );
      ia=[]; ia=size(FVCOM1.('s'));
     if length(size(ia))<=2,
      FVCOM1.('t') =permute(FVCOM1.('t')  ,ord3d);
      FVCOM1.('s') =permute(FVCOM1.('s')  ,ord3d);
    end
    FVCOM1.el= FVCOM.el;
    FVCOM1.h = FVCOM.h ;
    FVCOM1.xc= FVCOM.xc;
    FVCOM1.yc= FVCOM.yc;
    FVCOM1.u = permute(FVCOM.u ,ord3d );
    FVCOM1.v = permute(FVCOM.v ,ord3d );
    FVCOM1.ua= permute(FVCOM.ua,ord2d );
    FVCOM1.va= permute(FVCOM.va,ord2d );
    FVCOM1.ww= permute(FVCOM.ww,ord3d );
    FVCOM1.z = (FVCOM.z );
    FVCOM1.r = (FVCOM.r );
 %FVCOM1.vorticity    = permute(FVCOM.vorticity     );
  FVCOM1.short_wave   = permute(FVCOM.short_wave   ,ord2d );
  FVCOM1.net_heat_flux= permute(FVCOM.net_heat_flux,ord2d );
  FVCOM1.uwind_speed  = permute(FVCOM.uwind_speed  ,ord2d );
  FVCOM1.vwind_speed  = permute(FVCOM.vwind_speed  ,ord2d );
  FVCOM1.precip       = permute(FVCOM.precip       ,ord2d );
  FVCOM1.evap         = permute(FVCOM.evap         ,ord2d );

          out_mat=0;
  else

    FVCOM1.Date  =[FVCOM1.Date  ;Date];
    FVCOM1.s =[FVCOM1.s  ;single(FVCOM.s )];
    FVCOM1.t =[FVCOM1.t  ;single(FVCOM.t )];
    FVCOM1.el=[FVCOM1.el ;single(FVCOM.el)];
    FVCOM1.u =[FVCOM1.u  ;single(FVCOM.u )];
    FVCOM1.v =[FVCOM1.v  ;single(FVCOM.v )];
    FVCOM1.ua=[FVCOM1.ua ;single(FVCOM.ua)];
    FVCOM1.va=[FVCOM1.va ;single(FVCOM.va)];
    FVCOM1.ww=[FVCOM1.ww ;single(FVCOM.ww)];
    FVCOM1.z =[FVCOM1.z  ;single(FVCOM.z )];
    FVCOM1.r =[FVCOM1.r  ;single(FVCOM.r )];

  %FVCOM1.vorticity    =[FVCOM1.vorticity      ;single(FVCOM.vorticity     )];
  FVCOM1.short_wave    =[FVCOM1.short_wave     ;single(permute(FVCOM.short_wave,ord2d     ))];
  FVCOM1.net_heat_flux =[FVCOM1.net_heat_flux  ;single(permute(FVCOM.net_heat_flux,ord2d  ))];
  FVCOM1.uwind_speed   =[FVCOM1.uwind_speed    ;single(permute(FVCOM.uwind_speed,ord2d    ))];
  FVCOM1.vwind_speed   =[FVCOM1.vwind_speed    ;single(permute(FVCOM.vwind_speed,ord2d    ))];
  FVCOM1.precip        =[FVCOM1.precip         ;single(permute(FVCOM.precip,ord2d         ))];
  FVCOM1.evap          =[FVCOM1.evap           ;single(permute(FVCOM.evap ,ord2d          ))];

end