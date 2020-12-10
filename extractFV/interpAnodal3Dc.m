function [FPT] = interpAnodal3Dc(varargin)
% INTERPANODAL3DC - Interpolate quantities which are given at element
% corners to points at arbitrary spatial locations within those elements
%
% Note that the actual location of the new corner node is not required, as
% all distance calculations have been done already by the time this
% function is entered. What is read in is the identity of the node in the
% new arrays (ip) and a structure containing (amongst other things) a 
% matrix of precalculated distances/weights based on the nodes of the 
% containing element (nodeWeights), and the hydrodynamic model mesh.
%
% Inputs:   theMesh     - The original hydrodynamic model mesh
%           Field       - The field to be interpolated, expected to be an
%                           array of values (elements X depths, single time) 
%                           relating to each element centroid in mainMesh
%           nodeWeights - A structure containing distances and
%                           respective weights relating to centroids of the
%                           containing element in the original mesh, and
%                           its neighbours (max 4 elements total). This
%                           also contains the index of the containing
%                           element for the new location.
%           ip          - The index of the new point in the new mesh/arrays.
%           sigloc      - The depth of the point in sigma coordinates
%                           ([0,-1]) relative to the bathymetry depth in
%                           mainMesh.
%           lvls        - The depth of sigma levels in the hydrodynamic
%
% Outputs:  FPT         - interpolated value of the field at this location
%
% Authors:
%   Thomas Adams, SAMS
%   Dmitry Aleynik, SAMS
% 
% Revision history:
%   30/10/2018 Initial version produced under KTP project
%
    global runProperties

    nargs = length(varargin); 
    
    %disp('interpAnodal3Dc');
    %runProperties = varargin{1};
    theMesh = varargin{1};
    Field = varargin{2};
    nodeWeights = varargin{3};
    ip=varargin{4};
    % The 3D case - including depth location
    if nargs == 6
        sigloc =varargin{5};
        lvls  =varargin{6};
    end
    
    % Check whether to use CVM file or use calculated list of neighbours/weights
    invDist=0;
    if ~isempty(nodeWeights)
        %disp('Using nodeWeights');
        x0c = nodeWeights.x0c(ip);
        y0c = nodeWeights.y0c(ip);
        triId = nodeWeights.hydroElemIndex(ip);
        weights=nodeWeights.weights(ip,:);
        %if ~isempty(CVM)
        if any(strcmpi(runProperties.useCVM,{'Y','1'}))
            %disp('Using CVM');
            aw0=theMesh.aw0; % corner weighting of nodes for each element - always 0.33 per corner
            awx=theMesh.awx;
            awy=theMesh.awy;
            
            % ---------------- Reorder columns -----------------
            % *** Should no longer be required - fixed by Dima in metric
            % calculation ***
            %nv=theMesh.tri(:,[1,3,2]);
            nv=theMesh.tri;
            
            n1  = nv(triId,1); % The node vertices relating to each element - same as theMesh.trinodes 
            n2  = nv(triId,2);
            n3  = nv(triId,3);
        else
            invDist = 1;
            n1 = theMesh.tri(triId,1);
            n2 = theMesh.tri(triId,2);
            n3 = theMesh.tri(triId,3);
        end
    else
        error('CVM and weights are empty; cannot calculate interpolation');             
    end
        
    %% 2D case
    if nargs==4

        if invDist==1
            FPT = (weights(1)*Field(n1)+weights(2)*Field(n2)+weights(3)*Field(n3))...
                /(weights(1)+weights(2)+weights(3));
        else
            % ! linear interpolation of Field
            F0 = aw0(triId,1)*Field(n1)+aw0(triId,2)*Field(n2)+aw0(triId,3)*Field(n3);
            Fx = awx(triId,1)*Field(n1)+awx(triId,2)*Field(n2)+awx(triId,3)*Field(n3);
            Fy = awy(triId,1)*Field(n1)+awy(triId,2)*Field(n2)+awy(triId,3)*Field(n3);
            FPT = F0 + Fx*x0c + Fy*y0c;
            % F0 + Fx*xDistToCentre + Fy*yDistToCentre
%             disp(['oldElem ' num2str(i) ' nv ' num2str(n1) ' ' num2str(n2) ' ' num2str(n3), ...
%                 ' Field ' num2str(Field(n1)) ' ' num2str(Field(n2)) ' ' num2str(Field(n3)), ...
%                 ' F0 ' num2str(F0) ' Fx ' num2str(Fx) ' Fy ' num2str(Fy) ' FPT ' num2str(FPT)])
            disp('')
        end
        
    %% 3D case
    % Surrounding Element IDs
    elseif nargs==6

        if all([n1,n2,n3]) == 0
            FPT = Field(triId);
        else
            kb =length(theMesh.sigloc);
            kbm=kb; kbm1=kbm-1; kbm2=kb-2;

            zz1=theMesh.zz1;  
            z1=theMesh.z1; 
            dz1=theMesh.dz1; 
            dzz1=theMesh.dzz1;
        
            % ---------------------------          
            % Derive coefficients for vertical interpolation between the 
            % model sigma layers by determining the layer in which the
            % point resides
            
            if(lvls == kb)      
            % The case that the number of levels is the same as 
            % the number of hydrodynamic model layer boundaries
                if(sigloc >= z1(triId,1))
                    % disp(['top 3'])
                    k1 = 1;
                    k2 = 1;
                    alpha = -1;
                elseif(sigloc > z1(triId,kb))... then %!intermediate 
                    % disp([' intermediate 2']);
                    for k=1:kbm1
                        if(sigloc  < z1(triId,k) && sigloc >= z1(triId,k+1) )
                            k1 = k;
                            k2 = k+1 ;
                            alpha = (z1(triId,k)-sigloc)/dz1(triId,k);
                        end
                    end
                else
                    % !bottom
                    % disp([' bottom ']);
                    k1 = kbm1;
                    k2 = kbm1;
                    alpha = -1 ;
                end
            %elseif(lvls == kbm1)
            else
            % Originally, the case that the number of levels is the same as 
            % the number of hydrodynamic model vertical points. zz1 contains 11
            % vertical values on each row, the first 10 of which are
            % midpoints of the model layers.
            % Extended this to allow different numbers of layers
                if(sigloc >= zz1(triId,1))
                % Top layer ----
                    k1 = 1;
                    k2 = 1;
                    alpha = -1; % Enforces selection of a single layer
%                     disp(['sigloc=' num2str(sigloc) ' zz1(i,1)=' num2str(zz1(i,1)) ...
%                         ' alpha=' num2str(alpha)])
                elseif(sigloc > zz1(triId,kbm1)) 
                % Intermediate layer ----
                    for k=1:kbm2
                        if(sigloc  < zz1(triId,k) && sigloc >= zz1(triId,k+1) )
                            k1 = k;
                            k2 = k+1;
                            alpha = (zz1(triId,k)-sigloc)/dzz1(triId,k);
%                             disp(['sigloc=' num2str(sigloc) ' zz1(i,k)=' num2str(zz1(i,k)) ...
%                                 ' zz1(i,k+1)=' num2str(zz1(i,k+1)) ' alpha=' num2str(alpha)])
                            break;
                        end
                    end 
                else
                % Bottom layer
                    k1 = kbm1;
                    k2 = kbm1;
                    alpha = -1;
%                     disp(['sigloc=' num2str(sigloc) ' zz1(i,k1)=' num2str(zz1(i,k1)) ...
%                         ' alpha=' num2str(alpha)])
                end
%                 disp([ 'interpAnodal3D: invalid number of levels ',... 
%                     ' (must be equal to either kb or kbm1): lvls=' num2str(lvls) ]);
%                 return;
            end
            
            %% Tom's added inverse distance option
            if invDist==1
                f_upper = (weights(1)*Field(n1,k1)+weights(2)*Field(n2,k1)+weights(3)*Field(n3,k1))...
                    /(weights(1)+weights(2)+weights(3));
                
                if (k1==k2)
                    FPT=f_upper;
                else
                    f_lower = (weights(1)*Field(n1,k2)+weights(2)*Field(n2,k2)+weights(3)*Field(n3,k2))...
                        /(weights(1)+weights(2)+weights(3));
                    FPT=(alpha*f_lower)+((1.0-alpha)*f_upper);
                end
            else
        
                % !linear interpolation of Field
                f0 = aw0(triId,1)*Field(n1,k1)+aw0(triId,2)*Field(n2,k1)+aw0(triId,3)*Field(n3,k1);
                fx = awx(triId,1)*Field(n1,k1)+awx(triId,2)*Field(n2,k1)+awx(triId,3)*Field(n3,k1);
                fy = awy(triId,1)*Field(n1,k1)+awy(triId,2)*Field(n2,k1)+awy(triId,3)*Field(n3,k1);
                f_upper = f0 + fx*x0c + fy*y0c;

                if(k1 == k2)
                    FPT = f_upper;
                else

                    f0 = aw0(triId,1)*Field(n1,k2)+aw0(triId,2)*Field(n2,k2)+aw0(triId,3)*Field(n3,k2);
                    fx = awx(triId,1)*Field(n1,k2)+awx(triId,2)*Field(n2,k2)+awx(triId,3)*Field(n3,k2);
                    fy = awy(triId,1)*Field(n1,k2)+awy(triId,2)*Field(n2,k2)+awy(triId,3)*Field(n3,k2);
                    f_lower = f0 + fx*x0c + fy*y0c;

                    FPT = (alpha)*f_lower + (1.0-alpha)*f_upper ;
                    %disp('');
                end
            end
        end 
    end
end
