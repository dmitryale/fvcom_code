function [interpOut] = extractFiniteVolumeMethod_nopar(interpOut,hydroNodeInfo,newNodeInfo,newTri,hydroElemIndex)
% EXTRACTFINITEVOLUMEMETHOD - Extract the desired variables from the
% hydrodynamic model output files, and save into a Matlab structure. Once
% the method has completed, the structure contains both the subset of the
% original mesh output, and the values reinterpolated on to new locations.
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
%           newTri      - Information on new mehs traingle coordinates etc.
%           hydroElemIndex  - A list of the hydrodynamic model element
%                               containing each new node
%
% Outputs:  interpOut   - The reinterpolated data structure
%
% Authors:
%   Thomas Adams, SAMS
%   Dmitry Aleynik, SAMS
%
% Revision history:
%   30/10/2018 Initial version produced under KTP10483 project, MarineHarvest (Scotland) Ltd

    tic
    now2=now; 
    disp([ 'extractFiniteVolumeMethod_nopar start : ' datestr(now,31) ', dt=',num2str((now-now2)*24*3600) ' sec ' ]);    

    global runProperties mainMesh;
    %   global runProperties interpOut ;%#da  poolSize                       
    if 1
    %   global hydroNodeInfo newNodeInfo newTri hydroElemIndex nt nz vars_fv firstVar lastVar    
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
        [neighbours,~,weights,x0c,y0c]=getNeighbourElemDistsWeights(mainMesh,mainMesh.nbe,newNodeInfo.x,hydroElemIndex,1);
    else 
        [neighbours,~,weights,x0c,y0c]=getNeighbourElemDistsWeights(mainMesh,[],newNodeInfo.x,hydroElemIndex,1);
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
                 
        %#!da: get the elevations 
        [interpOut]= getElevationsNopar(interpOut,hydroNodeInfo,newNodeInfo,nt,nz,vars_fv,firstVar,nodeWeights);

        %% The T, S variables are recorded at nodes by FVCOM,                vars_fv(2,3)
        % The velocity variables are recorded at element centroids by FVCOM:vars_fv(4,5,6)
        % The will be interpolated to the nodes for the DEPOMOD files (that is,
        % all variables will end up on the same grid locations).   

        for vi=firstVar+1:lastVar ... % firstVar+2 % el,t,s firstVar=1;lastVar=6; %el,t,s ...

            va=char(vars_fv(vi));
            a=interpOut.(va);     % 2D or 3D
            % This statement sets the number of values searched to be at most
            % the size of the subset mesh - so comment out.
            sa=size(a); %sainity check what to use for extraction :
            if sa(end)==length(interpOut.nodes), cu=' interpAnodal3Dc'; useIntrp='Nodal'; end
            if sa(end)==length(interpOut.elems), cu=' interpAzonal3Dc'; useIntrp='Zonal'; end
            disp(['** ' va '  ** N_extractPts ' num2str(npv) ', [' num2str(size((a))) ...
                '], N_non-0_vals ' num2str(length(find(a))) ',' cu] );

            if any(strcmp(variablesToExtract,char(vars_fv(vi))))
                % The second and third variables
                %sva=[mainMesh.Nverts,nz];
                % Loop over time
                for it=1:nt,...
                    if vi <= firstVar+2  % t,s,
                        val=zeros(str2num(runProperties.hydroDepthLayers),mainMesh.Nverts);
                        val(:,hydroNodeInfo.nod)=double((a(it,:,:)));
                    else
                        val=zeros(str2num(runProperties.hydroDepthLayers),mainMesh.Nelems);
                        val(:,hydroNodeInfo.elm)=double(squeeze(a(it,:,:))); %3D u,v; 'tri_ids are the nearest'
                    end
                    aa=permute(val,[2,1]);

                    %% \\ // \\
                    %    Loop over spatial & vertical locations
                    for ip=1:npv

                        % disp(['location ' num2str(ip) ' --------------']);
                        % The first variable
                        for iz=1:nz
                            sig_loc = interpOut.(['fi_sig_loc' vab(1) ] )(it,iz,ip);

                            %lev=mainMesh.siglev; % As passed over
                            lev=newNodeInfo.sigloc;
                            if sig_loc==0
                                sig_loc=lev(1);
                            end
                            if isnan(sig_loc)
                                sig_loc=lev(end);
                            end

                            if vi<=firstVar+2 % t,s,
                                [val_1]=interpAnodal3Dc(mainMesh,aa,nodeWeights,ip,sig_loc,lvls);
                            else      % u,v,w
                                [val_1]=interpAzonal3Dc(mainMesh,aa,neighbourWeights,ip,sig_loc,lvls);
                            end

                            interpOut.(['fi_' va])(it,iz,ip)=val_1;% interpOutP(np,iz); %.(['fi_' va]);
                        end% iz
                    end %ip
                    %% \\ // \\
                end % it

%                 disp('vert sums at t=1')
%                 tmp=squeeze(interpOut.(['fi_' va])(1,:,:));
%                 for loc=1:npv
%                     disp(num2str(sum(tmp(:,loc))));
%                 end
%                 disp('depth sums at t=1')
%                 tmp=squeeze(interpOut.(['fi_' va])(1,:,:));
%                 for dep=1:nz
%                     disp(num2str(sum(tmp(dep,:))));
%                 end
%                 disp('time index sums at depth=1')
%                 tmp=squeeze(interpOut.(['fi_' va])(:,1,:));
%                 for it=1:nt
%                     disp(num2str(sum(tmp(it,:))));
%                 end
            end
        end %vi
        
 
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
     
    toc; 
    disp([ 'extractFiniteVolumeMethod_nopar End: ' datestr(now,31) ', dt=',num2str((now-now2)*24*3600) ' sec ' ]);    
    %  a1=squeeze(interpOut.fi_sig_loce(:,:,1)); plot(a1)
    %             interpOut_NoP=interpOut;
    %  a1=squeeze(interpOut_NoP.fi_u(:,:,1)); plot(a1,'linewidth',1); hold on;
end
