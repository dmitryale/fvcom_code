function [interpOut] = extractFiniteVolumeMethod_par(interpOut,hydroNodeInfo,newNodeInfo,newTri,hydroElemIndex)
% EXTRACTFINITEVOLUMEMETHOD - Extract the desired variables from the
% hydrodynamic model output files, and save into a Matlab structure. Once
% the method has completed, the structure contains both the subset of the
% original mesh output, and the values reinterpolated on to new locations.
%
% split va(riables) among pooSize (4) wokers
% elapse time (with npv=560 & pooSize=4) is 829.844 sec & 262.91 sec at 8 cpu
% java.util.concurrent.LinkedBlockingQueue (Java method) calls takes 99.5percet of time  
% r = q.poll(1,java.util.concurrent.TimeUnit.SECONDS); inside remoteparfor  L277
%
% Inputs:   interpOut     - The structure that will contain the interpolated
%                           data
%           hydroNodeInfo - Node and elements IDs(including boundary and
%                           coast) from the existing mesh that are within 
%                           the new domain (elmu), crossing the edge (elm),
%                           and within a threshold of the edge (elma). Also
%                           Boundary, coast, wet/dry nodes.
%           newNodeInfo - Node information to which the hydrodynamic data
%                           fields will be interpolated
%           newTri      - Information on new mesh traingle coordinates etc.
%           hydroElemIndex  - A list of the hydrodynamic model element
%                               containing each new node
%           nt          - Number of time records in hydrodynamic model file
%           nz          - Number of depth layers in hydrodynamic model file
%           vars_fv     - The list of variables for extraction (characters)
%           firstVar    - the numeric index of the first variable to
%                           extract (from vars_fv)
%           lastVar     - the numeric index of the last variable to
%                           extract (from vars_fv)
%
% Outputs:  interpOut   - The reinterpolated data structure
%
% Authors:
%   Thomas Adams, SAMS
%   Dmitry Aleynik, SAMS
%
% Revision history:
%   30/10/2018 Initial version produced under KTP10483 project, MarineHarvest (Scotland) Ltd
 
    now3=now; 
    disp([ 'extractFiniteVolumeMethod_par (va) start : ' datestr(now,31) ', dt=',num2str((now-now3)*24*3600) ' sec ' ]);    

    global runProperties mainMesh poolSize ;
 %  global runProperties interpOut   ;%#da                        
    if 1 % for debuging :
      %global hydroNodeInfo newNodeInfo newTri hydroElemIndex nt nz vars_fv firstVar lastVar    
        nt=str2num(runProperties.hydroTimeRecords);  
        firstVar=1; 
        lastVar=6;
        vars_fv = {'el';'t';'s'; 'u';'v';'w'; 'ua';'va';}; 
        nva = length(vars_fv);
        if newNodeInfo.sigloc(end)==-1
            nz=length(newNodeInfo.sigloc)-1;
        else
            nz=length(newNodeInfo.sigloc);
        end     
    end

    extractX = newNodeInfo.x(:,1); 
    extractY = newNodeInfo.x(:,2);
    xc = mainMesh.xc;
    yc = mainMesh.yc;
    lvls = nz;
    % Copy the interpolated coordinates to the output data structure
    interpOut.x = newNodeInfo.x;
    interpOut.lon = newNodeInfo.lon;
    interpOut.lat = newNodeInfo.lat;
    %firstVar=1; lastVar=6;      % el,t,s, u,v,ww  variables in vars_fv array
    vab=char(vars_fv(firstVar));% first variable to extract (elevation) - retained for use by later variables
    nva = length(vars_fv);
    
    % Check which variables to extract
    variablesToExtract=strsplit(runProperties.variables,',');

    %% Get the index of the main mesh element containing each node in the new mesh.
    % Need to identify which element of the main mesh (identified 
    % by nodeInfoExiting.elm) contains the spatial location
    % [cx,cy] in the new mesh
    % This is now calculated already in "compareMeshes"
    %hydroElemIndex=zeros(length(cx),1);
    for ip=1:length(extractX)
        %hydroElemIndex(ip)=whichElement(cx(ip),cy(ip),hydroNodeInfo.elm);
        % Catch the case where the centroid of the new element is outside
        % the old mesh
        if hydroElemIndex(ip)==0
            % cx and cy are lat/lon, so nearestCentroid uses use values in
            % mainMesh.uvnode, distance reported doesn't make too much sense
            [nearestID,dist] = nearestCentroid(newTri.geog(ip,1),newTri.geog(ip,2),mainMesh.xc(hydroNodeInfo.elm),mainMesh.yc(hydroNodeInfo.elm),hydroNodeInfo.elm);
            hydroElemIndex(ip) = nearestID;
            %warning(['New element ' int2str(ip) ' outside hydrodynamic mesh; using nearest element = ' num2str(nearestID) '(dist = ' num2str(dist) ')']);
        end
        %hydroElemIndex(ip)=tri_ids(ip);
    end
    
    %% Calculate distances from each new interpolation point to containing
    % element centroid and the bounding element centroids.
    % Higher powers (last argument) in this function give smoother
    % transition in values across elements; lower values make values more
    % homoegeneous within elements
    % i) Find neighbour elements
    
    % ii) Find distances and weights
    if any(strcmpi(runProperties.useCVM,{'Y','1'}))
        [neighbours,~,weights,x0c,y0c]=getNeighbourElemDistsWeights(mainMesh,[],newNodeInfo.x,hydroElemIndex,1);
    else
        [neighbours,~,weights,x0c,y0c]=getNeighbourElemDistsWeights(mainMesh,mainMesh.nbe,newNodeInfo.x,hydroElemIndex,1);
    end
    neighbourWeights=[];
    neighbourWeights.hydroElemIndex=hydroElemIndex;
    neighbourWeights.neighbours=neighbours;
    neighbourWeights.weights=weights;
    neighbourWeights.x0c=x0c;
    neighbourWeights.y0c=y0c;
    
    %% Calculate distances and weights from each new interpolation point
    % to the corner nodes of its containing element
    [distances,weights] = getNodeDistsWeights(mainMesh,newNodeInfo.x,hydroElemIndex,1);
    nodeWeights=[];
    nodeWeights.hydroElemIndex=hydroElemIndex;
    nodeWeights.distances=distances;
    nodeWeights.weights=weights;
    nodeWeights.x0c=x0c; % Add the centroid distances to nodeWeights too, as they were used in both zonal and nodal calculations.
    nodeWeights.y0c=y0c;
            
     %% The first three variables (elevation, temperature, salinity) are
    % recorded on nodes in the unstructured grid. Extract these first and
    % interpolate to new node locations for DEPOMOD file.
    
    % Pre-allocate arrays
                          npv=length(extractX);
%      if strcmpi(DRV,'C'), npv=11; end
    %disp('Rapid test of array sizes (no interpolation)');
    % var 1 (el)
    va=char(vars_fv(1));
    interpOut.(['fi_h_' va(1)])=zeros(npv,1);
    interpOut.(['fi_el' va(1)])=zeros(npv,nt); % 1
    interpOut.(['fi_sig_loc' va(1)])=zeros(nt,nz,npv);
    % var 2 (t)
    va=char(vars_fv(2));
    interpOut.(['fi_' va])=zeros(nt,nz,npv); % 2,3
    % var 3 (s)
    va=char(vars_fv(3));
    interpOut.(['fi_' va])=zeros(nt,nz,npv);
    % var 4 (u)
    va=char(vars_fv(4));
    interpOut.(['fi_' va])=zeros(nt,nz,npv);
    % var 5 (v)
    va=char(vars_fv(5));
    interpOut.(['fi_' va])=zeros(nt,nz,npv);
    % var 6 (w)
    va=char(vars_fv(6));
    interpOut.(['fi_' va])=zeros(nt,nz,npv);    
    
      
    test=0;
    if (test==0)
                 
        %#!da: get the el,
        [interpOut]= getElevationsNopar(interpOut,hydroNodeInfo,newNodeInfo,nt,nz,vars_fv,firstVar,nodeWeights);
        %  a1=squeeze(interpOut.fi_sig_loce(:,:,1)); plot(a1)
        nvars=lastVar-firstVar; 
        MxP = ceil(nvars/poolSize);  % MxP is the max number of calls to invoke parallel works 560/7=>80
        if MxP < 1, disp(['number of parralel calls < par.poolSize (' num2str(poolSize) '), replace run.usePar=0, return now !']); return; end

        %% The T, S variables are recorded at nodes by FVCOM,                vars_fv(2,3)
        % The velocity variables are recorded at element centroids by FVCOM:vars_fv(4,5,6)
        % The will be interpolated to the nodes for the DEPOMOD files (that is,
        % all variables will end up on the same grid locations).  
        BC = struct();
        for vi=firstVar+1:lastVar ...   
            va=char(vars_fv(vi));
            %       a=interpOut.(va);     % 2D or 3D
            BC(vi,1).va =va;                                     %-> the input to parloop
            BC(vi,1).vars(1:length(va)) =va;                     %-> the input to parloop
        end  

        MiP=zeros(poolSize,MxP);
        % we need re-alocate variables to match the number of wokers in a poolSize
        vi=0; ij=0;
        for jp=1:MxP  
            for kp=1:poolSize 
                vi=(firstVar+0+kp)+(jp-1)*poolSize;
                if vi > lastVar, continue ;end
                ij=ij+1;
                MiP(kp,jp)   = vi;   
                BC(vi,1).Mip = vi;
            end
        end
        % -------------
        % \\ // \\  
        for jp=1:MxP ,... % calls to invoke parfor kp cycles
            gm=MiP(:,jp);
            gm=gm(gm>0); % we will use only complete loops 
            Mip=gm;
            %        BC(vi,1).(['fi_sig_loc' vab(1) ] )(1:nt,1:nz,1:npv)=interpOut.(['fi_sig_loc' vab(1) ] )(1:nt,1:nz,1:npv);           
            % emptify local in/out structrure 
            interpOutP=[]; 
            for np=1:poolSize%+1
                interpOutP(np,1).out=zeros(nt,nz,npv); 
            end
            %% \\ // \\ // \\ // \\
            % parallel Loop itself: for each varible va we will use its own core-worker 
            % to go over [ temporal, vertical & spatial locations]      
            parfor np=1:poolSize+0 
                lp=np;%-1; % lp is the local counter at given processor np                     
                if lp > length(Mip), continue; end % skip empty va cycles
                vi = Mip(lp) ;
                if vi > lastVar || vi <=0 , continue; end % skip empty va cycles
                va = char(vars_fv(vi));
                a = interpOut.(va);                 % 2D or 3D
                b = interpOut.(['fi_sig_loc' vab(1) ]); 
                % This statement sets the number of values searched to be at most
                % the size of the subset mesh - so comment out.
                sa=size(a); %sainity check what to use for extraction : 
%                 if sa(end)==length(interpOut.nodes), cu=' interpAnodal3Dc'; useIntrp='Nodal'; end
%                 if sa(end)==length(interpOut.elems), cu=' interpAzonal3Dc'; useIntrp='Zonal'; end
                disp(['** ' va '  ** N_extractPts ' num2str(npv) ', [' num2str(size((a))) ...
                    '], N_non-0_vals ' num2str(length(find(a))) '; at woker # ' num2str(np) ', jp=' num2str(jp) ]);         
                %lev=mainMesh.siglev; % As passed over
                lev=newNodeInfo.sigloc;

                % clearvars jt mtvi go_fi vel vel_1 vel_2 vel_3 factm_a factm_b hh el1 el2 el3
                % Loop over space
                for ip=1:npv
                    % disp(['location ' num2str(ip) ' --------------']);
                    % Loop over time
                    for it=1:nt
                        for iz=1:nz                     
                            %sig_loc = interpOut.(['fi_sig_loc' vab(1) ] )(it,iz,ip);
                            sig_loc = b(it, iz,ip);

                            %lev=mainMesh.siglev; % As passed over
                            %lev=newNodeInfo.sigloc;
                            if sig_loc==0
                                sig_loc=lev(1);
                            end
                            if isnan(sig_loc)   
                                sig_loc=lev(end); 
                            end

                            if vi<=firstVar+2 % t,s,
                                val=zeros(str2num(runProperties.hydroDepthLayers),mainMesh.Nverts);
                                val(:,hydroNodeInfo.nod)=double(squeeze(a(it,:,:))); 
                                aa=permute(val,[2,1]);      
                                [val_1]=interpAnodal3Dc(runProperties,mainMesh,aa,nodeWeights,ip,sig_loc,lvls); 
                            else      % u,v,w
                                val=zeros(str2num(runProperties.hydroDepthLayers),mainMesh.Nelems);
                                val(:,hydroNodeInfo.elm)=double(squeeze(a(it,:,:))); %3D u,v; 'tri_ids are the nearest'
                                aa=permute(val,[2,1]);
                                [val_1]=interpAzonal3Dc(runProperties,mainMesh,aa,neighbourWeights,ip,sig_loc,lvls);                  
                            end

                            interpOutP(np,1).out(it,iz,ip)=val_1;       
                        end %iz
                    end %it
                end %ip
                interpOutP(np,1).va=va;  interpOutP(np,1).vi=vi;
            end % np %par-loop end
            %% \\ // \\ // \\ // \\
            %   collect all local np points from all wokers into the united global vi-array:              
            for np=1:poolSize+0   
                lp=np;%-1; % lp is the local counter at given processor np                     
                if lp > length(Mip), continue; end % skip empty va cycles
                vi = Mip(lp) ;
                if vi > lastVar || vi <=0 , continue; end % skip empty va cycles
                va = char(vars_fv(vi));          
                interpOut.(['fi_' va])=interpOutP(np,1).out;
%             for ip=1:npv                     
%                for it=1:nt     
%                 for iz=1:nz              
%           interpOut.(['fi_' va])(it,iz,ip)=interpOutP(np,1).out(it,iz,ip);                 
%                 end
%                end
%              end
            end                 
          
         end %jp MxP
        %% \\ // \\          
                         
     elseif test==1
        disp('Rapid test of array sizes (NO INTERPOLATION CARRIED OUT)');
        % var 1 (el)
        va=char(vars_fv(1));
        interpOut.(['fi_h_' va(1)])=ones(npv,1);
        interpOut.(['fi_el' va(1)])=ones(npv,nt); % 1
        interpOut.(['fi_sig_loc' va(1)])=ones(nt,nz,npv);
        % var 2 (t)
        va=char(vars_fv(2));
        interpOut.(['fi_' va])=ones(nt,nz,npv); % 2,3
        % var 3 (s)
        va=char(vars_fv(3));
        interpOut.(['fi_' va])=ones(nt,nz,npv);
        % var 4 (u)
        va=char(vars_fv(4));
        interpOut.(['fi_' va])=ones(nt,nz,npv);
        % var 5 (v)
        va=char(vars_fv(5));
        interpOut.(['fi_' va])=ones(nt,nz,npv);
        % var 6 (w)
        va=char(vars_fv(6));
        interpOut.(['fi_' va])=ones(nt,nz,npv);
        
    end % if
     
    disp([ 'extractFiniteVolumeMethod_par (va) End: ' datestr(now,31) ', dt=',num2str((now-now3)*24*3600) ' sec ' ]);    
% debuging test shows exact match with _nopar case: 
%               interpOut_Par=interpOut;
%  a2=squeeze(interpOut_Par.fi_u(:,:,1)); plot(a2,'linewidth',2); hold on;
%  a1=squeeze(interpOut_NoP.fi_u(:,:,1)); plot(a1,'-.','linewidth',1); hold on;

end
