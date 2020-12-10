function[Z,Z1,ZZ,ZZ1,DZ,DZ1,DZZ1,SIGMA] = setup_sigma_gen(fileToRead1,zl,Mesh )
%setup_sigma.m is based on fvcom's mod_setup.F subroutine sigma_geometric
% %INDEX_VERCOR, P_SIGMA, KB, DU2,DL2, N, M, NV, D);
% (c) dmitry.aleynik@sams.ac.uk 2011.03.28;                               &
%  __o_O__¬                                                               &
%  \_____/ ~~~~~~<@})))< ~ ~ ~~~~~ ~~~ SAMS                               &
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

   if ~exist('zl','var'),...
      zl=[0. 1.0  5.0 10.0 20.0 30.0      50.0     100.      150. 200.  225.];  % 11
   end

  if ~exist('fileToRead1','var')
%  fileToRead1='sigma_gen.dat'; % general
   fileToRead1='sigma.dat';     % tanh
%  fileToRead1='sigma_u.dat';   % uniform
  end

%  deafult value :
   P_SIGMA=1.0 ;

  fid=fopen([fileToRead1]);
   a=[]; a=fgetl(fid); clear b k;
  if length(a)>=5 && strcmp(a(1:6),'NUMBER'), k=find(a=='='); b=a(k+1:end); KB=str2double(b);end
  clear a b k  TYPE ;
      j=0;
  for i=1:19,...
     a=[]; a=fgetl(fid);
     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=5 && strcmp(a(1:5),'SIGMA'), k=find(a=='='); b=a(k+1:end); TYPE=strtrim(b); j=1 ; end
     clear b k ;
%% =================== 0
   if j==1 && length(TYPE) >=7 && strcmp(TYPE(1:7),'UNIFORM') ==1 ,...

   end %j=1
%% =================== I
   if j==1 && length(TYPE) >=9 && strcmp(TYPE(1:9),'GEOMETRIC')==1 ,...
     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=12 && strcmp(a(1:12),'SIGMA POWER'), k=find(a=='='); b=a(k+1:end);  P_SIGMA=str2double(b); end
     clear b k ;          %  deafult value :
   end %j=1
%% =================== II
   if j==1 && length(TYPE) >= 4 && strcmp(TYPE(1:4),'TANH')==1 ,...
     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=6 && strcmp(a(1:2),'DU'), k=find(a=='='); b=a(k+1:end); DU=str2double(b); end
     clear b k ;
     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=6 && strcmp(a(1:2),'DL'), k=find(a=='='); b=a(k+1:end); DL=str2double(b); end
     clear b k ;
   end %j=1
%% =================== III
   if j==1 && length(TYPE) >= 11 && strcmp(TYPE(1:11),'GENERALIZED') == 1 ,...
     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=6 && strcmp(a(1:2),'DU'), k=find(a=='='); b=a(k+1:end); DU=str2double(b); end
     clear b k ;
     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=6 && strcmp(a(1:2),'DL'), k=find(a=='='); b=a(k+1:end); DL=str2double(b); end
     clear b k ;
     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=6 && strcmp(a(1:2),'MI'), k=find(a=='='); b=a(k+1:end); HMIN1=str2double(b); end % MinConstantDepth'
     clear b k ;
     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=6 && strcmp(a(1:2),'KU'), k=find(a=='='); b=a(k+1:end); KU=str2double(b); end
    clear b k ;
     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=6 && strcmp(a(1:2),'KL'), k=find(a=='='); b=a(k+1:end); KL=str2double(b); end
    clear b k ;
     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=6 && strcmp(a(1:3),'ZKU'), k=find(a=='='); b=a(k+1:end); ZKU=str2num(b); end

     if isempty(a)==0 && a(1:1)~='!' && ...
        length(a)>=6 && strcmp(a(1:3),'ZKL'), k=find(a=='='); b=a(k+1:end); ZKL=str2num(b); end
   end %j

  end
  clear a b i k;
 %================
  SIGMA.KB = KB;  clear KB;% # OF SIGMA LEVELS(KB)
  SIGMA.TYPE=TYPE;

  if length(TYPE)>=7 && strcmp(TYPE(1:7),'UNIFORM') ==1,...
  INDEX_VECTOR =1 ;
  SIGMA.P_SIGMA = P_SIGMA  ; %=1
  end

  if length(TYPE)>=9 && strcmp(TYPE(1:9),'GEOMETRIC') ==1, ...
  SIGMA.P_SIGMA = P_SIGMA  ; %=2
  INDEX_VECTOR =1 ;
  end

  if length(TYPE)>=4 && strcmp(TYPE(1:4),'TANH')==1,...
  SIGMA.DU  =DU  ; %DU2
  SIGMA.DL  =DL  ; %DL2
  clear DU DL;
  INDEX_VECTOR =2 ;
  end

 if length(TYPE)>=11 && strcmp(TYPE(1:11),'GENERALIZED')==1,...
  SIGMA.DU  =DU  ; %DUU
  SIGMA.DL  =DL  ; %DLL
  SIGMA.HMIN1 =HMIN1 ;
  SIGMA.KU  =KU  ;
  SIGMA.KL  =KL  ;
  SIGMA.ZKU =ZKU ;
  SIGMA.ZKL =ZKL ;
  clear DU DL MCD KU KL ZKU ZKL HMIN1;
  INDEX_VECTOR =3 ;
 end

  SIGMA.IDX=INDEX_VECTOR   ;
  SIGMA.P_SIGMA = P_SIGMA  ;

fclose(fid);
%
%v2.7.1
%  INDEX_VERCOR=2; % INDEX_VERCOR=1;%   INDEX_VERCOR=2; % INDEX_VERCOR=3;
%     P_SIGMA=2.0 ;  %P_SIGMA=1.0 ;
 % DU2 = 0.001; DL2 = 0.001; %DEFAULT equal; %DU2 = 0.200; DL2 = 2.000; % condesed near bottom
 % DU2 = 2.000; DL2 = 2.000; % condesed near surface and bottom

   DPTHSL = zl;          % m
      KSL = length(zl);   % number of Z standard levels
SIGMA.KSL = KSL;
SIGMA.DPTHSL = DPTHSL;
      KB   = SIGMA.KB;    % number of sigma levels, must be odd KB=15,11


  N =Mesh.Nelems; % number of elements (triangles' centroids)
  M =Mesh.Nverts; % number of nodes    (vertices <| )
  NV=Mesh.tri   ;  % NV(1:N,3);   % NODE NUMBERING FOR ELEMENTS (0:NT)
  D =Mesh.x(:,3);    % depth at nodes  noe elevation ;
  H =Mesh.x(:,3);    % depth at nodes


 KBM1= KB-1;
   ZTMP=zeros(KB,1);
      Z=zeros(M,KB);
     Z1=zeros(N,KB);
     D3n=zeros(M,KB);
     D3c=zeros(N,KB);
     Dc =zeros(N,1);

IDX = SIGMA.IDX ;% INDEX_VECTOR;

%% case (1):   ===============
if IDX==1, ...

     P_SIGMA=SIGMA.P_SIGMA;

  if P_SIGMA == 1,...
   for K=1:KB, ...
       ZTMP(K,1) = -((K-1)/(KB-1))^P_SIGMA ;
   end
  else
   for K=1:(KB+1)/2,...
       ZTMP(K) = -((K-1)/((KB+1)/2-1) )^P_SIGMA/2 ;
   end
   for K=((KB+1)/2)+1 : KB,...
       ZTMP(K) =  ((KB-K)/((KB+1)/2-1) )^P_SIGMA/2 - 1.0 ;
   end
  end

   for J=1:M,...
     for K=1:KB,...
       Z(J,K)=ZTMP(K);
     D3n(J,K)=D(J)*Z(J,K);  %depth at nods   DEPTH_Z(K)=Z(I,K)*D(I) !LAYER CENTER DEPTH
     end
   end

   for I=1:N,...
       Dc(I)  =(D(NV(I,1) )+D(NV(I,2)  )+ D(NV(I,3)  )) / 3.0 ;
     for K=1:KB,...
       Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+ Z(NV(I,3),K)) / 3.0 ;
      D3c(I,K)=Dc(I)*Z1(I,K);  %depth at centroids
     end
   end
end

%%  case (ii) Tanh
if IDX==2,...
  DL2 =SIGMA.DL;
  DU2 =SIGMA.DU;

   for K=1:KBM1,...
     X1=DL2+DU2;
     X1=X1*(KBM1-K)/KBM1;
     X1=X1-DL2;
     X1=tanh(X1);
     X2=tanh(DL2);
     X3=X2+tanh(DU2);
     for J=1:M,...
        Z(J,K+1)=(X1+X2)/X3 - 1.0;
     end
     for I=1:N,...
       Z1(I,K+1)=(X1+X2)/X3 - 1.0;
     end
   end

   for K=1:KB,...
     for J=1:M,...
     D3n(J,K)=D(J)*Z(J,K);  %depth at nods !LAYER CENTER DEPTH
     end

      for I=1:N,...
    if K==1,    Dc(I)=(D(NV(I,1) )+D(NV(I,2)  )+ D(NV(I,3)  )) / 3.0 ;  end
      D3c(I,K)=Dc(I)*Z1(I,K);  %depth at centroids
      end
   end

end
%% =============III
if IDX==3,...
  DUU=    SIGMA.DU  ;
  DLL=    SIGMA.DL  ;
 HMIN1 = SIGMA.HMIN1 ;
 KU = SIGMA.KU  ;
 KL = SIGMA.KL  ;
 ZKU= SIGMA.ZKU ;
 ZKL= SIGMA.ZKL ;

    for I=1:M,...
         if(H(I) < HMIN1),...
            Z(I,1)=0.0;

            DL2=0.001;
            DU2=0.001;

            for K=1:KBM1, ...
               X1=DL2+DU2;
               X1=X1*(KBM1-K)/KBM1;
               X1=X1-DL2;
               X1=tanh(X1);
               X2=tanh(DL2);
               X3=X2+tanh(DU2);

               Z(I,K+1)=(X1+X2)/X3-1.0;
            end
        else
            DR=(H(I)-DUU-DLL)/H(I)/(KB-KU-KL-1);

            Z(I,1)=0.0;

            for K=2:KU+1,...
               Z(I,K)=Z(I,K-1)-ZKU(K-1)/H(I);
            end

            for K=KU+2:KB-KL,...
               Z(I,K)=Z(I,K-1)-DR;
            end

            KK=0;
            for K=KB-KL+1:KB,...
               KK=KK+1;
               Z(I,K)=Z(I,K-1)-ZKL(KK)/H(I);
            end
         end
    end

      for I=1:N,...
         for K=1:KB,...
            Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+Z(NV(I,3),K))/3.0 ;
         end
      end
end% ! SIGMA_GENERALIZED

%% ===================

%% derivatives for all three cases:
for K=1:KB-1,...
     for I=1:M,...
       ZZ(I,K)  = 0.5*(Z(I,K)+Z(I,K+1));
       DZ(I,K)  = Z(I,K)-Z(I,K+1);
     end
     for I=1:N,...
       DZ1(I,K)  = Z1(I,K)-Z1(I,K+1);
       ZZ1(I,K)  = 0.5*(Z1(I,K)+Z1(I,K+1));
     end
end
   for I=1:M,...
     ZZ(I,KB) = 2.0*ZZ(I,KB-1)-ZZ(I,KB-2);
   end
   for I=1:N,...
     ZZ1(I,KB) = 2.0*ZZ1(I,KB-1)-ZZ1(I,KB-2);
   end

   for K=1:KB-1 ,...
     for I=1:N,...
       DZZ1(I,K) = ZZ1(I,K)-ZZ1(I,K+1);
     end
   end
    DZZ1(:,KB) = 0.0;
end
%% --------------

% %==============================