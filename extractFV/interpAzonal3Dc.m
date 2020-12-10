function [interpVal] = interpAzonal3Dc(varargin)
% INTERPAZONAL3DC - Interpolate values from hydrodynamic mesh element
% centroids to corners nodes of new mesh.
%
% Note that the actual location of the new corner node is not required, as
% all distance calculations have been done already by the time this
% function is entered. What is read in is the identity of the node in the
% new arrays (ip) and a structure containing (amongst other things) a 
% matrix of precalculated distances/weights based on the centroids of 
% neighbouring elements (neighbourWeights), and (optionally, can be empty) 
% a structure containing general metrics relating to the original mesh 
% elements (CVM).
%
% Inputs:   mainMesh    - The original hydrodynamic model mesh
%           Field       - The field to be interpolated, expected to be an
%                           array of values (elements X depths, single time) 
%                           relating to each element centroid in mainMesh
%           CVM         - Continuous Volume Metrics structure describing
%                           features of the original mesh (optional)
%           neighbourWeights - A structure containing distances and
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
% Outputs:  interpVal   - The interpolated value at the specified location.
%
% Authors:
%   Thomas Adams, SAMS
%   Dmitry Aleynik, SAMS
% 
% Revision history:
%   30/10/2018 Initial version produced under KTP10483 project, MarineHarvest (Scotland) Ltd
%

    global runProperties

    nargs = length(varargin); 
    
    %runProperties = varargin{1};
    theMesh = varargin{1};
    Field = varargin{2};
    neighbourWeights = varargin{3};
    ip=varargin{4};
    % The 3D case - including depth location
    if nargs == 6
        sigloc =varargin{5};
        lvls  =varargin{6};
    end

    % Check whether to use CVM file or use calculated list of neighbours/weights
    invDist=0;
    if ~isempty(neighbourWeights)
        %disp('Using neighbourWeights');
        x0c = neighbourWeights.x0c(ip);
        y0c = neighbourWeights.y0c(ip);
        triId = neighbourWeights.hydroElemIndex(ip);
        weights=neighbourWeights.weights(ip,:);
        if any(strcmpi(runProperties.useCVM,{'Y','1'}))
            %disp('Using CVM');
            a1u=theMesh.a1u;
            a2u=theMesh.a2u;
            
            % -------------- Reorder columns -------------- 
            % Note here, altered ordering of last three indices to
            % work with Dmitry's indexing (possibly a result of
            % Fortran/Matlab indexing differences)
            % *** Should no longer be required - fixed by Dima in metric
            % calculation ***
            %nbe=theMesh.nbe(:,[1,3,2]); 
            nbe=theMesh.nbe;
                        
            e1  = nbe(triId,1);
            e2  = nbe(triId,2);
            e3  = nbe(triId,3);
        else
            invDist = 1;
            a1u=[];
            a2u=[];
            e1 = neighbourWeights.neighbours(ip,2);
            e2 = neighbourWeights.neighbours(ip,3);
            e3 = neighbourWeights.neighbours(ip,4);
        end
        disp('');
    else
        error('CVM and weights are empty; cannot calculate interpolation');             
    end
    %invDist

    %% 2D
    % Surrounding Element IDs
    if nargs==4

        if all([e1,e2,e3]) == 0
            interpVal = Field(triId);
        else
            if invDist==1
                interpVal = (weights(1)*Field(triId)+weights(2)*Field(e1)+weights(3)*Field(e2)+weights(4)*Field(e3))...
                    /(weights(1)+weights(2)+weights(3)+weights(4));                
            else  
                % interpolate 2D Field to the location
                Fx = a1u(triId,1)*Field(triId)+a1u(triId,2)*Field(e1)+a1u(triId,3)*Field(e2)+a1u(triId,4)*Field(e3);
                Fy = a2u(triId,1)*Field(triId)+a2u(triId,2)*Field(e1)+a2u(triId,3)*Field(e2)+a2u(triId,4)*Field(e3);
                interpVal = Field(triId) + (Fx*x0c) + (Fy*y0c);
            end
        end
        
    %% 3D
    %Surrounding Element IDs
    elseif nargs==6

        kb=length(theMesh.sigloc);
        kbm=kb; 
        kbm1=kbm-1;     
        kbm2=kb-2;
        
        zz1=theMesh.zz1;  
        z1=theMesh.z1; 
        dz1=theMesh.dz1; 
        dzz1= theMesh.dzz1;

        % ---------------------------     ip=1;    
        % ! Determine the layer in which the point resides
        if(lvls == kb) % then       
            % !top
            if(sigloc >= z1(triId,1))
                %disp(['top 3'])
                k1 = 1;
                k2 = 1;
                alpha = -1;
            elseif(sigloc > z1(triId,kb)) % then %!intermediate 
                %disp([' intermediate 2']);
                for k=1:kbm1,...
                    if(sigloc  < z1(triId,k) && sigloc >= z1(triId,k+1) )
                        k1 = k;
                        k2 = k+1 ;
                        alpha = (z1(triId,k)-sigloc)/dz1(triId,k);
                    end
                end
            else
                % !bottom
                %disp([' bottom ']);
                k1 = kbm1;
                k2 = kbm1;
                alpha = -1 ;
            end
        %elseif(lvls == kbm1)
        else
        % Altered to allow arbitrary number of levels - just find the
        % nearest layer centre
        % !top
            if(sigloc >= zz1(triId,1))
                k1 = 1;
                k2 = 1;
                alpha = -1;
                %disp(['top 1'])
            elseif(sigloc > zz1(triId,kbm1))  % !intermediate 
                for k=1:kbm2
                    if(sigloc  < zz1(triId,k) && sigloc >= zz1(triId,k+1) ), ...
                        k1 = k;
                        k2 = k+1 ;
                        alpha = (zz1(triId,k)-sigloc)/dzz1(triId,k);
                        break ;% return; exit
                    end
                    %disp([' intermediate 1'])
                end 
            else
                % ! top
                %disp(['top 2'])
                k1 = kbm1;
                k2 = kbm1;
                alpha = -1;
            end
%         else           
%             disp([ 'interp_azonal3D: invalid number of levels ',... &
%                  ' (must be equal to either kb or kbm1): lvls=' num2str(lvls) ]);
%             return;
        end
    %% ===================
        if all([e1,e2,e3]) == 0 % assign zero weight to coastal edge element :
            % Field(1,:)=0;
            if e1==0, e1=1; end
            if e2==0, e2=1; end
            if e3==0, e3=1; end
        end

        % Tom inverse distance weighting - to compare against
        % existing a1u, a2u scheme that produces strange results
         
        if invDist==1
            f_upper = (weights(1)*Field(triId,k1)+weights(2)*Field(e1,k1)+weights(3)*Field(e2,k1)+weights(4)*Field(e3,k1))...
                /(weights(1)+weights(2)+weights(3)+weights(4));

            if (k1==k2)
                interpVal=f_upper;
            else
                f_lower = (weights(1)*Field(triId,k2)+weights(2)*Field(e1,k2)+weights(3)*Field(e2,k2)+weights(4)*Field(e3,k2))...
                    /(weights(1)+weights(2)+weights(3)+weights(4));
                interpVal=(alpha*f_lower)+((1.0-alpha)*f_upper);
            end
        else  
            % interpolate Field to the location based on values in CVM.a1u, CVM.a2u
            fx=(a1u(triId,1)*Field(triId,k1))+(a1u(triId,2)*Field(e1,k1))+(a1u(triId,3)*Field(e2,k1))+(a1u(triId,4)*Field(e3,k1));
            fy=(a2u(triId,1)*Field(triId,k1))+(a2u(triId,2)*Field(e1,k1))+(a2u(triId,3)*Field(e2,k1))+(a2u(triId,4)*Field(e3,k1));
            f_upper=Field(triId,k1)+(fx*x0c)+(fy*y0c);

            if (k1==k2)
                interpVal=f_upper;
            else
                fx=a1u(triId,1)*Field(triId,k2)+a1u(triId,2)*Field(e1,k2)+a1u(triId,3)*Field(e2,k2)+a1u(triId,4)*Field(e3,k2);
                fy=a2u(triId,1)*Field(triId,k2)+a2u(triId,2)*Field(e1,k2)+a2u(triId,3)*Field(e2,k2)+a2u(triId,4)*Field(e3,k2);
                f_lower=Field(triId,k2)+(fx*x0c)+(fy*y0c);

                interpVal=(alpha*f_lower)+((1.0-alpha)*f_upper);
            end 
        end
        disp('');
    end
end