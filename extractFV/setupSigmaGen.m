function [Z,Z1,ZZ,ZZ1,DZ,DZ1,DZZ1,SIGMA] = setupSigmaGen()
% SETUPSIGMAGEN Generate sigma layer depths across entire mesh, based on a
% defined distribution of depth levels (paramters contained within a file
% to be read in). These are added to the global mesh variable, so the
% output values are not required under normal usage. The function is based 
% on fvcom's mod_setup.F subroutine sigma_geometric.
%
% In FVCOM output, velocities are reported at element centroids, and
% scalars (e.g. temperature, salinity) are reported at corder nodes.
%
% Inputs:   []
%
% Outputs:  Z           - Node sigma layer depths(proportion)
%           Z1          - Centroid sigma layer depths (proportion)
%           ZZ          - Node sigma layer midpoint depths (proportion)
%           ZZ1         - Centroid sigma layer midpoint depths (proportion)
%           DZ          - Node sigma layer thickness (proportion)
%           DZ1         - Centroid sigma layer thickness (proportion)
%           DZZ1        - Centroid midpoint differences (proportion)
%           SIGMA       - Structure containing details of Sigma layer
%                           structure parameters
%
% Also created  D3n     - sigma layer depths at all element corner nodes (in metres)
%               D3c     - sigma layer depth at all element centroids (in metres)
%               Dc      - depth at element centroids (in metres)

    disp('Setting up input depth layers');

    global runProperties mainMesh;
    
    N =mainMesh.Nelems; % number of elements (triangles' centroids)
    M =mainMesh.Nverts; % number of nodes    (vertices <| )
    NV=mainMesh.trinodes;  % NV(1:N,3);   % NODE NUMBERING FOR ELEMENTS (0:NT)
    D =mainMesh.x(:,3);    % depth at nodes  noe elevation ;
    H =mainMesh.x(:,3); 
    KB=str2num(runProperties.hydroDepthLayers)+1;
    
    % A file potentially containing information on sigma-layers; note that
    % this may not exist
    fileToRead=[runProperties.hydroMeshDir '/sigma.dat'];
    
    % First check whether the mesh already contains the relevant
    % information; SSM meshes should
    if isfield(mainMesh,'siglay') && isfield(mainMesh,'siglev') 
        % Check that this has the correct format and size
        if(size(mainMesh.siglay,2)==1)
            mainMesh.siglay=repmat(mainMesh.siglay',mainMesh.Nverts,1);
            mainMesh.siglev=repmat(mainMesh.siglev',mainMesh.Nverts,1);
        end
        
        if ~isequal(size(mainMesh.siglay),[mainMesh.Nverts,str2num(runProperties.hydroDepthLayers)+1])        
            error('Sigma layer boundaries defined in mesh file but not correct size (Nverts X runProperties.hydroDepthLayers)');
        end
        if ~isequal(size(mainMesh.siglev),[mainMesh.Nverts,str2num(runProperties.hydroDepthLayers)+1])        
            error('Sigma layer midpoints defined in mesh file but not correct size (Nverts X runProperties.hydroDepthLayers)');
        end
        Z=mainMesh.siglev;
        Z1=repmat(mainMesh.siglev(1,:),mainMesh.Nelems,1);
        ZZ=mainMesh.siglay;
        ZZ1=repmat(mainMesh.siglev(1,:),mainMesh.Nelems,1);
        %DZ=
        %DZ1=
        %DZZ1=
            
    elseif exist(fileToRead, 'file') == 2
        % Do stuff based on reading the file in here

        % Set the number of sigma levels (should be first line of file)
        fid=fopen([fileToRead]);
        a=[]; 
        a=fgetl(fid); 
        clear b k;
        if length(a)>=5 && strcmp(a(1:6),'NUMBER')
            k=find(a=='='); 
            b=a(k+1:end); 
            KB=str2double(b);
        end
        clear a b k  TYPE ;
        j=0;

        P_SIGMA=1.0;
        % Read subsequent lines in turn to get sigma layer parameters.
        % TODO: Clean up this code to remove looping and improve readability
        for i=3:19,...
            a=[]; 
            a=fgetl(fid); 
            if isempty(a)==0 && a(1:1)~='!' && length(a)>=5 && strcmp(a(1:5),'SIGMA')
                k=find(a=='='); 
                b=a(k+1:end); 
                TYPE=strtrim(b); 
                j=1; 
            end;
            clear b k ;
            %% =================== 0          
            if j==1 && length(TYPE) >=7 && strcmp(TYPE(1:7),'UNIFORM') ==1 ,...           

            end %j=1      
            %% =================== I          
            if j==1 && length(TYPE) >=9 && strcmp(TYPE(1:9),'GEOMETRIC')==1 ,...           
                if isempty(a)==0 && a(1:1)~='!' && length(a)>=12 && strcmp(a(1:12),'SIGMA POWER') 
                    k=find(a=='='); 
                    b=a(k+1:end);  
                    P_SIGMA=str2double(b); 
                end
                clear b k ;          %  deafult value :
            end %j=1  
            %% =================== II          
            if j==1 && length(TYPE) >= 4 && strcmp(TYPE(1:4),'TANH')==1 ,...           
                if isempty(a)==0 && a(1:1)~='!' && length(a)>=6 && strcmp(a(1:2),'DU') 
                    k=find(a=='='); 
                    b=a(k+1:end); 
                    DU=str2double(b); 
                end;
                clear b k ;     
                if isempty(a)==0 && a(1:1)~='!' && length(a)>=6 && strcmp(a(1:2),'DL')
                    k=find(a=='='); 
                    b=a(k+1:end);
                    DL=str2double(b); 
                end;
                clear b k ;
            end %j=1  
            %% =================== III     
            if j==1 && length(TYPE) >= 11 && strcmp(TYPE(1:11),'GENERALIZED') == 1 ,...           
                if isempty(a)==0 && a(1:1)~='!' && length(a)>=6 && strcmp(a(1:2),'DU'), 
                    k=find(a=='='); 
                    b=a(k+1:end); 
                    DU=str2double(b); 
                end
                clear b k ;
                if isempty(a)==0 && a(1:1)~='!' && length(a)>=6 && strcmp(a(1:2),'DL')
                    k=find(a=='=');
                    b=a(k+1:end); 
                    DL=str2double(b); 
                end
                clear b k ;
                if isempty(a)==0 && a(1:1)~='!' && length(a)>=6 && strcmp(a(1:2),'MI')
                    k=find(a=='='); 
                    b=a(k+1:end); 
                    HMIN1=str2double(b); 
                end % MinConstantDepth'
                clear b k ;
                if isempty(a)==0 && a(1:1)~='!' && length(a)>=6 && strcmp(a(1:2),'KU')
                    k=find(a=='=');
                    b=a(k+1:end); 
                    KU=str2double(b); 
                end
                clear b k ;
                if isempty(a)==0 && a(1:1)~='!' && length(a)>=6 && strcmp(a(1:2),'KL')
                    k=find(a=='='); 
                    b=a(k+1:end); 
                    KL=str2double(b); 
                end
                clear b k ;
                if isempty(a)==0 && a(1:1)~='!' && length(a)>=6 && strcmp(a(1:3),'ZKU')
                    k=find(a=='=');
                    b=a(k+1:end); 
                    ZKU=str2num(b); 
                end

                if isempty(a)==0 && a(1:1)~='!' && length(a)>=6 && strcmp(a(1:3),'ZKL')
                    k=find(a=='='); 
                    b=a(k+1:end); 
                    ZKL=str2num(b); 
                end
            end %j

        end
        clear a b i k;
        %================
        SIGMA.KB = KB;  % # OF SIGMA LEVELS(KB)
        SIGMA.KSL = KB; % # OF Z-LEVELS(KB)
        clear KB;
        SIGMA.TYPE=TYPE; 

        if length(TYPE)>=7 && strcmp(TYPE(1:7),'UNIFORM') ==1,...  
            INDEX_VECTOR =1 ;  
            SIGMA.P_SIGMA = 1.0  ; %=1
        end;   

        if length(TYPE)>=9 && strcmp(TYPE(1:9),'GEOMETRIC') ==1, ...  
            SIGMA.P_SIGMA = P_SIGMA  ; %=2
            INDEX_VECTOR =1 ;         
        end;   

        if length(TYPE)>=4 && strcmp(TYPE(1:4),'TANH')==1,...  
            SIGMA.DU  =DU  ; %DU2
            SIGMA.DL  =DL  ; %DL2
            clear DU DL;
            INDEX_VECTOR =2 ;  
        end;   

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
        end;

        SIGMA.IDX=INDEX_VECTOR   ;
        SIGMA.P_SIGMA = P_SIGMA  ; 

        fclose(fid);
        %
        %v2.7.1
        %  INDEX_VERCOR=2; % INDEX_VERCOR=1;%   INDEX_VERCOR=2; % INDEX_VERCOR=3;
        %     P_SIGMA=2.0 ;  %P_SIGMA=1.0 ; 
         % DU2 = 0.001; DL2 = 0.001; %DEFAULT equal; %DU2 = 0.200; DL2 = 2.000; % condesed near bottom
         % DU2 = 2.000; DL2 = 2.000; % condensed near surface and bottom   

        % Commented out these three lines 11/06/18 as dont actually use the z
        % levels read in - just the number. This is contained within the file
        % sigma.dat.
        %DPTHSL = zl;          % m
        %5KSL = length(zl);   % number of Z standard levels 
        %SIGMA.DPTHSL = DPTHSL; 

        KB   = SIGMA.KB;    % number of sigma levels, must be odd KB=15,11

        % depth at nodes  


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
                end;
            else 
                for K=1:(KB+1)/2,...
                    ZTMP(K) = -((K-1)/((KB+1)/2-1) )^P_SIGMA/2 ;
                end
                for K=((KB+1)/2)+1 : KB,...
                    ZTMP(K) =  ((KB-K)/((KB+1)/2-1) )^P_SIGMA/2 - 1.0 ;
                end
            end;  

            for J=1:M,...
                for K=1:KB,...
                    Z(J,K)=ZTMP(K);
                    D3n(J,K)=D(J)*Z(J,K);  %depth at nodes   DEPTH_Z(K)=Z(I,K)*D(I) !LAYER CENTER DEPTH
                end;
            end;

            for I=1:N,...
                Dc(I)  =(D(NV(I,1) )+D(NV(I,2)  )+ D(NV(I,3)  )) / 3.0 ; 
                for K=1:KB,...
                    Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+ Z(NV(I,3),K)) / 3.0 ;  
                    D3c(I,K)=Dc(I)*Z1(I,K);  %depth at centroids
                end;
            end;   
        end;   

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
                end;               
                for I=1:N,...
                    Z1(I,K+1)=(X1+X2)/X3 - 1.0;
                end;
            end;

            for K=1:KB,...
                for J=1:M,...
                    D3n(J,K)=D(J)*Z(J,K);  %depth at nods !LAYER CENTER DEPTH  
                end;     

                for I=1:N,...
                    if K==1,    
                        Dc(I)=(D(NV(I,1) )+D(NV(I,2)  )+ D(NV(I,3)  )) / 3.0 ;  
                    end;
                    D3c(I,K)=Dc(I)*Z1(I,K);  %depth at centroids
                end;
            end;  

        end;
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
                    end;

                    for K=KU+2:KB-KL,...
                       Z(I,K)=Z(I,K-1)-DR;
                    end

                    KK=0;
                    for K=KB-KL+1:KB,...
                       KK=KK+1;
                       Z(I,K)=Z(I,K-1)-ZKL(KK)/H(I);
                    end               
                 end
            end;

              for I=1:N,...
                 for K=1:KB,...
                    Z1(I,K)=(Z(NV(I,1),K)+Z(NV(I,2),K)+Z(NV(I,3),K))/3.0 ;
                 end
              end
        end;% ! SIGMA_GENERALIZED

    end

    %% ===================

    %% derivatives for all three cases:
    for K=1:KB-1,...
        % Nodes
        for I=1:M,...
            % Layer midpoints
            ZZ(I,K)  = 0.5*(Z(I,K)+Z(I,K+1));
            % Layer thickness
            DZ(I,K)  = Z(I,K)-Z(I,K+1);
        end;
        % Centroids
        for I=1:N,...  
            % Layer midpoints
            ZZ1(I,K)  = 0.5*(Z1(I,K)+Z1(I,K+1));
            % Layer thickness
            DZ1(I,K)  = Z1(I,K)-Z1(I,K+1);
        end;
    end;
    % Node bottom layer midpoint
    for I=1:M,...  
        ZZ(I,KB) = 2.0*ZZ(I,KB-1)-ZZ(I,KB-2);
    end;
    % Centroid bottom layer midpoint
    for I=1:N,...  
        ZZ1(I,KB) = 2.0*ZZ1(I,KB-1)-ZZ1(I,KB-2);
    end;  
   
    % Difference in midpoints
    for K=1:KB-1 ,...
        for I=1:N,...
            DZZ1(I,K) = ZZ1(I,K)-ZZ1(I,K+1);
        end
    end
    DZZ1(:,KB) = 0.0;
    
    %% Add calculated variables to mainMesh
    mainMesh.z   = Z;
    mainMesh.zz  = ZZ;
    mainMesh.z1  = Z1;          
    mainMesh.zz1 = ZZ1; 
    mainMesh.dz1 = DZ1;
    mainMesh.dzz1= DZZ1; 
    
    mainMesh.sigloc=Z1(1,:)'; % Sigma layer depth boundaries
    mainMesh.siglev=Z1(1,:)'; % Sigma layer depth boundaries
    mainMesh.siglay=ZZ1(1,:)'; % Sigma layer depth midpoints
    
end
