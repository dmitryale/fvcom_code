function [FPT] = interp_anodal3Dc(varargin)
%(c)dmitry.aleynik@sams.ac.uk, 2011.08.22 , based on mod_utils.F fvcom
% [FPT] = interp_anodal(xloc,yloc,            tri_id, Field,Mesh, CVM); % 2D        &
% [FPT] = interp_anodal(xloc,yloc,sigloc,lvls,tri_id, Field); % 3D        &
% Edited : Dmitry Aleynik , 23/03/2017 %!da offlag.F subroutine interp_kh &
%  __o_O__¬                                                               &
%  \_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS , AZIMUTH 2011-2013           &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

     nargs = length(varargin) ; narg = nargs-4;
%    disp ( [' narg=', num2str(narg) ])

 % global Mesh CVM
  Mesh = varargin{narg+1};
  CVM  = varargin{narg+2};

aw0=CVM.aw0;
awx=CVM.awx;
awy=CVM.awy;
 nv=CVM.NV;
%   xc=Mesh.uvnode; %(:,1);
%   yc=xc(:,2);  xc=xc(:,1);
%   yc=Mesh.uvnode(:,2);
  xc = varargin{narg+3};
  yc = varargin{narg+4};

          xloc   = varargin{1};
          yloc   = varargin{2};
  if narg==4 , i = varargin{3}; end;
  if narg==6 , i = varargin{5}; end;

 % !offset from element center:  in meters
 if abs(yloc)<90 ,... % defined(SPHERICAL)
     rearth    = 6371.0e03 ;      % !! earth radius in meters
     deg2rad   = pi*2/360.0;      % !! radians/degree
     tpi       = deg2rad*rearth;  % !tpi=pi*rearth/180.=3.14159265/180.0*6371.*1000. =111.19km

     y0c = tpi * (yloc - yc( i)) ; %2397m
         dx_sph =  xloc - xc(i) ;
     if (dx_sph > 180.0 ),
         dx_sph = dx_sph -360.0 ;
     elseif (dx_sph < -180.0 ),
        dx_sph = dx_sph + 360.0 ;
     end;
     x0c = tpi * cos(deg2rad*(yloc + yc(i))*0.5 ) * dx_sph ; %-1160 m
 else
     x0c = xloc - xc(i);
     y0c = yloc - yc(i);
 end;

 if narg==4 ,...
             tri_id = varargin{3};
              Field = varargin{4};
  i=tri_id;
% Fieled = t, s ,el, h
%    !element location (i) and surrounding nodes (n1,n2,n3)
     n1  = nv(i,1);
     n2  = nv(i,2);
     n3  = nv(i,3);

       if all([n1,n2,n3]) == 0
        FPT = Field(i);
       else
      %if all([n1,n2,n3]) ~=0 ,...
%      FPT=(Field(n1) + Field(n2) + Field(n3) ) /3;
% x0c = xloc - xc(i);
% y0c = yloc - yc(i);

% ! linear interpolation of Field
     F0 = aw0(i,1)*Field(n1)+aw0(i,2)*Field(n2)+aw0(i,3)*Field(n3);
     Fx = awx(i,1)*Field(n1)+awx(i,2)*Field(n2)+awx(i,3)*Field(n3);
     Fy = awy(i,1)*Field(n1)+awy(i,2)*Field(n2)+awy(i,3)*Field(n3);
     FPT = F0 + Fx*x0c + Fy*y0c;
       end
%
elseif narg==6
%% 3D case
% Surrounding Element IDs
%   disp(['narg=',num2str(narg)] );

%     xloc  =varargin{1};
%     yloc  =varargin{2};
    sigloc  =varargin{3};
    lvls  =varargin{4};
    tri_id=varargin{5};
    Field =varargin{6};
%     Field=Field';
          i=tri_id;

     n1  = nv(i,1);
     n2  = nv(i,2);
     n3  = nv(i,3);

  if all([n1,n2,n3]) == 0,
        FPT = Field(i); % FPT=Field(n1,k1);
   else

%          N=Mesh.Nelems; % w=79244;
%          M=Mesh.Nverts; % ww=46878;
%        get_siglocs; % setup_sigma_gen;
          kb =length(Mesh.sigloc); %! da correct 11
          kbm=kb; kbm1=kbm-1; kbm2=kb-2; % 10
%               z1=repmat(Mesh.sigloc,[1,M ])';  % nodes 11
%              zz1=repmat(Mesh.siglay,[1,N ])';  % elms 10
        zz1=Mesh.zz1;  z1=Mesh.z1; dz1=Mesh.dz1; dzz1= Mesh.dzz1;

% ---------------------------     ip=1;

%      ! Determine the layer in which the point resides
     if(lvls == kbm1),....
%         !top
        if(sigloc >= zz1(i,1)),
           k1 = 1;
           k2 = 1;
           alpha = -1;
%                      disp(['top 1'])
        elseif(sigloc > zz1(i,kbm1)),  % !intermediate
        for k=1:kbm2,
              if(sigloc  < zz1(i,k) && sigloc >= zz1(i,k+1) ), ...
                 k1 = k;
                 k2 = k+1 ;
                 alpha = (zz1(i,k)-sigloc)/dzz1(i,k);
                 break ;% return; exit
              end;
%           disp([' intermediate 1'])
        end

        else
%            ! top
%                     disp(['top 2'])
           k1 = kbm1;
           k2 = kbm1;
           alpha = -1;
        end;
      elseif(lvls == kb),... then
%         !top
       if(sigloc >= z1(i,1)),... then
%                         disp(['top 3'])
           k1 = 1;
           k2 = 1;
           alpha = -1;
        elseif(sigloc > z1(i,kb)),... then %!intermediate
%                               disp([' intermediate 2']);
      for k=1:kbm1,...
              if(sigloc  < z1(i,k) && sigloc >= z1(i,k+1) ),... then
                 k1 = k;
                 k2 = k+1 ;
                 alpha = (z1(i,k)-sigloc)/dz1(i,k);
              end;
      end;
        else
%            !bottom
 %     disp([' bottom ']);
           k1 = kbm1;
           k2 = kbm1;
           alpha = -1 ;
        end;
     else

       disp([ 'interp_anodal3D: invalid number of levels ',... &
             ' (must be equal to either kb or kbm1): lvls=' num2str(lvls) ]);
         return;
     end;

%      !linear interpolation of Field
     f0 = aw0(i,1)*Field(n1,k1)+aw0(i,2)*Field(n2,k1)+aw0(i,3)*Field(n3,k1);
     fx = awx(i,1)*Field(n1,k1)+awx(i,2)*Field(n2,k1)+awx(i,3)*Field(n3,k1);
     fy = awy(i,1)*Field(n1,k1)+awy(i,2)*Field(n2,k1)+awy(i,3)*Field(n3,k1);
     f_upper = f0 + fx*x0c + fy*y0c;

     if(k1 == k2),... then
        FPT = f_upper;

     else

        f0 = aw0(i,1)*Field(n1,k2)+aw0(i,2)*Field(n2,k2)+aw0(i,3)*Field(n3,k2);
        fx = awx(i,1)*Field(n1,k2)+awx(i,2)*Field(n2,k2)+awx(i,3)*Field(n3,k2);
        fy = awy(i,1)*Field(n1,k2)+awy(i,2)*Field(n2,k2)+awy(i,3)*Field(n3,k2);
        f_lower = f0 + fx*x0c + fy*y0c;

        FPT = (alpha)*f_lower + (1.0-alpha)*f_upper ;
     end;

  end;
 end;
end
