function [output] = interp_azonal3Dc(varargin)
% [output] = interp_azonal3Dv2(xloc,yloc,tri_id,Field);             %2D   &
% [output] = interp_azonal3Dv2(xloc,yloc,sigloc,lvls,tri_id,Field); %3D   &
% (c)dmitry.aleynik@sams.ac.uk, 2011.08.22;% based on FVCOM part.tracking &
% Updated: Meghan Rochford, 09/01/2017;KTP10483,MarineHarvest(Scotland)Ltd&
% Edited : Dmitry Aleynik , 21/04/2017  %!da                 mod_utils.F  &
%  __o_O__¬                                                               &
%  \_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS , AZIMUTH 2011-2013           &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

     nargs = length(varargin) ; narg = nargs-4;
%    disp ( [' narg=', num2str(narg) ])

% global Mesh CVM
% load('Mesh.mat'); % load('CVM_090101.mat');
  Mesh = varargin{narg+1};
% CVM  = varargin{narg+2};

a1u= Mesh.a1u;
a2u= Mesh.a2u;
nbe= Mesh.NBE;
  xc = varargin{narg+3};
  yc = varargin{narg+4};
%   xc=Mesh.uvnode; %(:,1);
%   yc=xc(:,2);
%   xc=xc(:,1);
              xloc=varargin{1};
              yloc=varargin{2};
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

%% 2D
% Surrounding Element IDs
if narg==4
    % disp(['narg=',num2str(narg)] );
    tri_id=varargin{3};
    Field=varargin{4};

    i=tri_id;

    e1  = nbe(i,1);
    e2  = nbe(i,2);
    e3  = nbe(i,3);

    if all([e1,e2,e3]) == 0
        output = Field(i);
    else
%         x0c = xloc - xc(i);
%         y0c = yloc - yc(i);

        % interpolate 2D Field to the location
        Fx = a1u(i,1)*Field(i)+a1u(i,2)*Field(e1)+a1u(i,3)*Field(e2)+a1u(i,4)*Field(e3);
        Fy = a2u(i,1)*Field(i)+a2u(i,2)*Field(e1)+a2u(i,3)*Field(e2)+a2u(i,4)*Field(e3);
        output = Field(i) + (Fx*x0c) + (Fy*y0c);
    end

%% 3D
%Surrounding Element IDs
elseif narg==6

%      disp(['narg=',num2str(narg)] );
%    xloc =varargin{1};
%    yloc =varargin{2};
   sigloc =varargin{3};
    lvls  =varargin{4};
    tri_id=varargin{5};
    Field =varargin{6};
   % Field =Field' ;
          i=tri_id;

    e1  = nbe(i,1);
    e2  = nbe(i,2);
    e3  = nbe(i,3);


%          N=Mesh.Nelems; % w=79244;
%          M=Mesh.Nverts; % ww=46878;
%        get_siglocs; % setup_sigma_gen;
         kb=length(Mesh.sigloc); %! da correct 11
         kbm=kb; kbm1=kbm-1;     kbm2=kb-2; % 10
%               z1=repmat(Mesh.sigloc,[1,M ])'; % nodes 11
%              zz1=repmat(Mesh.siglay,[1,N ])'; % elms 10
        zz1=Mesh.zz1;  z1=Mesh.z1; dz1=Mesh.dz1; dzz1= Mesh.dzz1;

% ---------------------------     ip=1;

%      ! Determine the layer in which the point resides
     if(lvls == kbm1),....
%         !top
        if(sigloc >= zz1(i,1)),
           k1 = 1;
           k2 = 1;
           alpha = -1;
%                       disp(['top 1'])
        elseif(sigloc > zz1(i,kbm1)),  % !intermediate
        for k=1:kbm2,
              if(sigloc  < zz1(i,k) && sigloc >= zz1(i,k+1) ), ...
                 k1 = k;
                 k2 = k+1 ;
                 alpha = (zz1(i,k)-sigloc)/dzz1(i,k);
                 break ;% return; exit
              end;
%            disp([' intermediate 1'])
        end

        else
%            ! top
%                      disp(['top 2'])
           k1 = kbm1;
           k2 = kbm1;
           alpha = -1;
        end;
      elseif(lvls == kb),... then
%          !top
       if(sigloc >= z1(i,1)),... then
%                          disp(['top 3'])
           k1 = 1;
           k2 = 1;
           alpha = -1;
        elseif(sigloc > z1(i,kb)),... then %!intermediate
%                                  disp([' intermediate 2']);
       for k=1:kbm1,...
              if(sigloc  < z1(i,k) && sigloc >= z1(i,k+1) ),... then
                 k1 = k;
                 k2 = k+1 ;
                 alpha = (z1(i,k)-sigloc)/dz1(i,k);
              end;
       end;
        else
%            !bottom
%       disp([' bottom ']);
           k1 = kbm1;
           k2 = kbm1;
           alpha = -1 ;
        end;
     else

       disp([ 'interp_azonal3D: invalid number of levels ',... &
             ' (must be equal to either kb or kbm1): lvls=' num2str(lvls) ]);
         return;
     end;
%% ===================
       if all([e1,e2,e3]) == 0,... % assigne zero weight to coastal edge element :
                     % Field(1,:)=0;
           if e1==0, e1=1; end
           if e2==0, e2=1; end
           if e3==0, e3=1; end
       end;

                % interpolate Field to the location
  fx=(a1u(i,1)*Field(i,k1))+(a1u(i,2)*Field(e1,k1))+(a1u(i,3)*Field(e2,k1))+(a1u(i,4)*Field(e3,k1));
  fy=(a2u(i,1)*Field(i,k1))+(a2u(i,2)*Field(e1,k1))+(a2u(i,3)*Field(e2,k1))+(a2u(i,4)*Field(e3,k1));
  f_upper=Field(i,k1)+(fx*x0c)+(fy*y0c);

                if (k1==k2)
                    output=f_upper;

                else
  fx=a1u(i,1)*Field(i,k2)+a1u(i,2)*Field(e1,k2)+a1u(i,3)*Field(e2,k2)+a1u(i,4)*Field(e3,k2);
  fy=a2u(i,1)*Field(i,k2)+a2u(i,2)*Field(e1,k2)+a2u(i,3)*Field(e2,k2)+a2u(i,4)*Field(e3,k2);
  f_lower=Field(i,k2)+(fx*x0c)+(fy*y0c);

                    output=(alpha*f_lower)+((1.0-alpha)*f_upper);
                end

%       end;

end;
end

