function[interpOut ]= getElevationsNopar(interpOut,hydroNodeInfo,newNodeInfo,nt,nz,vars_fv,firstVar,nodeWeights)
% GETELEVATIONSNOPAR - Extract elevations from the hydrodynamic output file
% and interpolate to all sigma layers in the output fileds for consumption
% by NewDEPOMOD.
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
%           nt          - Number of time records in hydrodynamic model file
%           nz          - Number of depth layers in hydrodynamic model file
%           vars_fv     - The list of variables for extraction (characters)
%           firstVar    - the numeric index of the first variable to
%                           extract (from vars_fv)
%           lastVar     - the numeric index of the last variable to
%                           extract (from vars_fv)
%
% Outputs:  interpOut        - interpolated value of the field at this location
%
% Authors:
%   Thomas Adams, SAMS
%   Dmitry Aleynik, SAMS
% 
% Revision history:
%   30/10/2018 Initial version produced under KTP project

    global mainMesh  

    extractX = newNodeInfo.x(:,1); 
    npv=length(extractX);
 
    % The first variable iz el :     
    %   for vi=firstVar:firstVar % el,t,s firstVar=1;lastVar=6; %el,t,s ...
    vi=firstVar;
    va=char(vars_fv(vi));           
    a=interpOut.(va);                  
    sa=size(a);        %sainity check what to use for extraction : 
    if sa(end)==length(interpOut.nodes), cu=' interpAnodal3Dc'; end
    if sa(end)==length(interpOut.elems), cu=' interpAzonal3Dc'; end
         
    disp(['** ' va ' ** N_extractPts ' num2str(npv) ', dimensions are [' num2str(size((a))) ...
        ']; N_non-0_vals ' num2str(length(find(a))), cu ]);
    % Loop over spatial locations   
    for ip=1:npv                 
        % Loop over time
        for it=1:nt
            %fprintf('-- ip %d hydroElemIndex(ip) %d vi %d it %d --\n', ip, hydroElemIndex(ip), vi, it);
            if it==1 % get the depths
                % Case that no bathymetry is provided: interpolate
                % existing bathymetry values 
                [hh] =interpAnodal3Dc(mainMesh,mainMesh.depth,nodeWeights,ip); %2D

                if ~isfield(newNodeInfo,'depth')
                    interpOut.(['fi_h_' va(1)])(ip,1)=hh;
                    % Case that new mesh contains bathymetry: use these
                    % values as provided
                else
                    interpOut.(['fi_h_' va(1)])(ip,1) = newNodeInfo.depth(ip);
                end

            else
                hh=interpOut.(['fi_h_' va(1)])(ip,1) ;
            end

            elev = zeros(1,mainMesh.Nverts);
            elev(1,hydroNodeInfo.nod) = interpOut.el(it,: );

            [el1]=interpAnodal3Dc(mainMesh,squeeze(elev(1,:)),nodeWeights,ip);

            interpOut.(['fi_el' va(1)])(ip,it)=el1 ;
            %z_tmp = -mainMesh.siglev*(hh+el1);
            z_tmp = newNodeInfo.sigloc*(hh+el1);

            % Work out the location of depth levels
            for iz=1:nz ...
                % Use the calculated depth at output sigma layer to query the
                % depths of the hydrodynamic model mesh at the same
                % location, and obtain a hydrodynamic model sigma
                % layer value. 
                % TODO: Could switch the vertical
                % interpolation part of the main interp routines to
                % live here, create a structure first time around
                % and use for other 5 variables.

                %disp('----------- debugging values ------------------')    
                %iz
                %z_tmp(iz)
                lev=mainMesh.siglev*hh;
                %mainMesh.siglev
                sig_loc=interp1(lev,mainMesh.siglev,z_tmp(iz));

                if sig_loc==0
                    sig_loc=lev(1);   
                end
                if isnan(sig_loc)
                    %sig_loc=lev(end);
                    sig_loc=-1;
                end

                interpOut.(['fi_sig_loc' va(1)])(it,iz,ip)=sig_loc;
            end %iz
        end%it
    end% ip
end