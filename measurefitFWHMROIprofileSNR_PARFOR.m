function     [dist, Int, fitG,ContourSegment,meshSegment, distO,c, CentreFit, xi, yi, IntInt, subsegmentArea, subsegmentVol, radii, RatioR, messedupframes,NoiseIdx] = measurefitFWHMROIprofileSNR_PARFOR(OriginalImage,pos,subNormals,mesh,meshSpacing, ContnextXsnake, ContnextYsnake,extent,Elength, StartFrame, EndFrame, segmentedMitoIDImageOld, MidPoint)

messedupframes=[];
Image = OriginalImage;
  % Image{f}= im2uint16(Image{f});
     
    I(:,:,1)=Image;
    %Icrop(:,:,f)=imcrop(I(:,:,f),(pos./3)); 
    
    subNormals(all(subNormals==0,2), :)=[];  
    TotDist=sqrt((subNormals(:,1)-subNormals(:,3)).*(subNormals(:,1)-subNormals(:,3))+(subNormals(:,2)-subNormals(:,4)).*(subNormals(:,2)-subNormals(:,4)));
    
%     subNormals=cellfun(@(x) x*px,subNormals1,'un',0);

%     p
% waitbar((i-StartFrame)/(EndFrame-StartFrame));
    if Elength<size(extent,1)
        if extent(Elength,1)>size(subNormals,1) 
            ContourSegment=[subNormals(end:extent(Elength,2),1) subNormals(end:extent(Elength,2),2) subNormals(end:extent(Elength,2),3) subNormals(end:extent(Elength,2),4)];
            meshSegment=[mesh(end:extent(Elength,2),1) mesh(end:extent(Elength,2),2) mesh(end:extent(Elength,2),3) mesh(end:extent(Elength,2),4)];
            MidPoints=[MidPoint(end:extent(Elength,2),1) MidPoint(end:extent(Elength,2),2)];
        elseif extent(Elength,2)>size(subNormals,1)
            ContourSegment=[subNormals(extent(Elength,1):end,1) subNormals(extent(Elength,1):end,2) subNormals(extent(Elength,1):end,3) subNormals(extent(Elength,1):end,4)];
            meshSegment=[mesh(extent(Elength,1):end,1) mesh(extent(Elength,1):end,2) mesh(extent(Elength,1):end,3) mesh(extent(Elength,1):end,4)];
            MidPoints=[MidPoint(extent(Elength,1):end,1) MidPoint(extent(Elength,1):end,2)];
        else
            ContourSegment=[subNormals(extent(Elength,1):extent(Elength,2),1) subNormals(extent(Elength,1):extent(Elength,2),2) subNormals(extent(Elength,1):extent(Elength,2),3) subNormals(extent(Elength,1):extent(Elength,2),4)];
            meshSegment=[mesh(extent(Elength,1):extent(Elength,2),1) mesh(extent(Elength,1):extent(Elength,2),2) mesh(extent(Elength,1):extent(Elength,2),3) mesh(extent(Elength,1):extent(Elength,2),4)];
            MidPoints=[MidPoint(extent(Elength,1):extent(Elength,2),1) MidPoint(extent(Elength,1):extent(Elength,2),2)];
        end
    else
        ContourSegment=[subNormals(extent(end,1):extent(end,2),1) subNormals(extent(end,1):extent(end,2),2) subNormals(extent(end,1):extent(end,2),3) subNormals(extent(end,1):extent(end,2),4)];
        meshSegment=[mesh(extent(end,1):extent(end,2),1) mesh(extent(end,1):extent(end,2),2) mesh(extent(end,1):extent(end,2),3) mesh(extent(end,1):extent(end,2),4)];
        MidPoints=[MidPoint(extent(end,1):extent(end,2),1) MidPoint(extent(end,1):extent(end,2),2)]
    end
    
    for j=1:length(ContourSegment)
        [cx{j,1}, cy{j,1}, c{j,1}, xi{j,1}, yi{j,1}]=improfile(I(:,:,1), [ContourSegment(j,1), ContourSegment(j,3)], [ContourSegment(j,2),ContourSegment(j,4)], 'bilinear' );
    end
% end
% close(h);
 
[SignalIdx, NoiseIdx]=measurefitFWHM_SNR(1, 1, segmentedMitoIDImageOld, cx, cy, ContourSegment);


% h = waitbar(0,'Finding valleys and peaks...');

% for i=StartFrame:EndFrame
%     waitbar((i-StartFrame)/(EndFrame-StartFrame));

    for j=1:length(ContourSegment)
        Int{j,1}=c{j,1};
        if min(Int{j,1})<0
%         MinI{j,1}=min(c{j,1});
            Int{j,1}=c{j,1}-min(c{j,1});
        else
            Int{j,1}=c{j,1}; %background subtracted profile
        end
        Int{j,1}(isnan(Int{j,1}(:,1)),:)=[];
        % Int{j,1}=lssmooth(Int{j,1},4);

        %transform real space coordinates to distance along line and compute
        %integrated intensity for each profile
        Dx{j,1}=cx{j,1}-xi{j,1}(1,1);
        Dy{j,1}=cy{j,1}-yi{j,1}(1,1);
        distO{j,1}= sqrt((Dx{j,1}.*Dx{j,1})+(Dy{j,1}.*Dy{j,1})); %this is the x distance along x-axis
        dist{j,1}= sqrt((Dx{j,1}.*Dx{j,1})+(Dy{j,1}.*Dy{j,1})); %this is the x distance along x-axis
        dist{j,1}(isnan(dist{j,1}(:,1)),:)=[];
        
%         if size(Int{j,1},1)==size(dist{j,1},1);
            IntInt(j,:)=trapz(dist{j,1}, Int{j,1}); %integrated intensity for each profile
%         end

        %calculate truncated cone area
        radii=sqrt((meshSegment(:,1)-meshSegment(:,3)).^2+(meshSegment(:,2)-meshSegment(:,4)).^2)/2;
        
        for r=1:(length(meshSegment)-1)
            subsegmentArea(r,:)=(pi()*(radii(r,1)+radii(r+1,1))*sqrt((radii(r,1)-radii(r+1,1))^2+(meshSpacing)^2));
            subsegmentVol(r, :)=(pi()*(1/3)*(meshSpacing)*((radii(r,1))^2+(radii(r+1,1))^2)+radii(r+1,1)*radii(r,1));
            RatioR(r, :)=radii(r+1,1)/radii(r,1);
        end

        % find peaks 
        [pks{j,1},locs{j,1},width{j,1}, prominence{j,1}] = findpeaks(Int{j,1},dist{j,1}); %pks and locs are the Y and X values of the peaks
        [~,IDXpks{j,1},~]=intersect(dist{j,1},locs{j,1});

            % find peaks in SIGNAL
            IDXpksS{j,1}=intersect(IDXpks{j,1},SignalIdx{j,1});

            % find peaks in NOISE
            IDXpksN{j,1}=intersect(IDXpks{j,1},NoiseIdx{j,1});

        % find valleys
        [valpks{j,1},vallocs{j,1},valwidth{j,1}, valprominence{j,1}] = findpeaks(-Int{j,1},dist{j,1}); %pks and locs are the Y and X values of the peaks
        [~, IDXval{j,1}, ~]=intersect(dist{j,1},vallocs{j,1});
        
            % find valleys in SIGNAL
            IDXvalS{j,1}=intersect(IDXval{j,1},SignalIdx{j,1});

            % find valleys in NOISE
            IDXvalN{j,1}=intersect(IDXval{j,1},NoiseIdx{j,1});
            
       % find CENTRAL peak
       
            % if signal peaks exist
            if isempty(IDXpksS{j,1})==0
                % find main peak %MidPoint
                d1=sqrt((ContourSegment(j,3)-MidPoints(j,1))^2 + (ContourSegment(j,4)-MidPoints(j,2))^2);
                d=sqrt((ContourSegment(j,1)-ContourSegment(j,3))^2 + (ContourSegment(j,2)-ContourSegment(j,4))^2);
                pkdist{j,1}=(d1/d)*max(dist{j,1});
                knnsearch(locs{j,1},pkdist{j,1},'k',1);
                Maxrow(:,j)=knnsearch(locs{j,1},pkdist{j,1},'k',1);
                PkInt{j,1}=pks{j,1}(Maxrow(:,j));
                Pkdist{j,1}=locs{j,1}(Maxrow(:,j));
%                 PkID{j,1}=knnsearch(dist{j,1},pkdist{j,1},'k',1);

                PkID{j,1}=find(dist{j,1}==locs{j,1}(Maxrow(:,j)));

%                 PkInt{j,1}=???
%                 Pkdist{j,1}=???
%                 PkID{j,1}=???
                
            % find valley height boundaries
                if isempty(IDXvalS{j,1})==0
                    VHup(j,1)=max(Int{j,1}(IDXvalS{j,1}))/PkInt{j,1};
                else
                    VHup(j,1)=1;
                end
                if isempty(IDXvalN{j,1})==0
%                     Int{j,1}(IDXvalN{j,1})
%                     max(Int{j,1}(IDXvalN{j,1}))
%                     PkInt{j,1}
                    VHlow(j,1)=max(Int{j,1}(IDXvalN{j,1}))/PkInt{j,1};
                else
                    VHlow(j,1)=0;
                end
                if VHup(j,1)>VHlow(j,1)
                    VH(j,1)=(VHup(j,1)+VHlow(j,1))/2;
                else
                    VH(j,1)=0.99*min([VHup(j,1) VHlow(j,1)]);
                end
            
            % find peak height boundaries
                if isempty(IDXpksS{j,1})==0
                    PHlow(j,1)=max(Int{j,1}(IDXpksS{j,1}))/PkInt{j,1};
                else
                    PHlow(j,1)=1;
                end
%                 if isempty(IDXpksN{j,1})==0
%                     PHup(j,1)=max(Int{j,1}(IDXpksN{j,1}))/PkInt{j,1};
%                 else
                    PHup(j,1)=max(Int{j,1}(NoiseIdx{j,1}))/PkInt{j,1};
%                 end
                if PHup(j,1)>PHlow(j,1)
                    PH(j,1)=(PHup(j,1)+PHlow(j,1))/2;
                else
                    PH(j,1)=1.01*max([PHup(j,1) PHlow(j,1)]);
                end
            else
                PkInt{j,1}=max(Int{j,1});
                Pkdist{j,1}=mean(dist{j,1});
                PkID{j,1}=round(mean(1:length(dist{j,1})));
                PH(j,1)=1.01*max(Int{j,1})/PkInt{j,1};
                VH(j,1)=0;
            end
%             [VH(j,1) PH(j,1)]
        % cut away the noise valleys
            
        IDXvalR{j,1}=IDXvalN{j,1}(dist{j,1}(IDXvalN{j,1})<dist{j,1}(PkID{j,1}));
        IDXvalL{j,1}=IDXvalN{j,1}(dist{j,1}(IDXvalN{j,1})>dist{j,1}(PkID{j,1}));
        
        IDXvalR{j,1}=IDXvalR{j,1}(Int{j,1}(IDXvalR{j,1})<VH(j,1)*PkInt{j,1});
        IDXvalL{j,1}=IDXvalL{j,1}(Int{j,1}(IDXvalL{j,1})<VH(j,1)*PkInt{j,1});
        IDXvalR{j,1}=max(IDXvalR{j,1});
        IDXvalL{j,1}=min(IDXvalL{j,1});
        
        if isempty(IDXvalR{j,1})==1
            IDXvalR{j,1}=1;
        end
        if isempty(IDXvalR{j,1})==1
            IDXvalL{j,1}=length(dist{j,1});
        end
        
        if IDXvalR{j,1}<IDXvalL{j,1}
            IDXfit{j,1}=[IDXvalR{j,1}:IDXvalL{j,1}];
        else
            IDXfit{j,1}=SignalIdx{j,1};
        end
        
        % cut away the noise peaks
                
        IDXpksR{j,1}=IDXpksN{j,1}(dist{j,1}(IDXpksN{j,1})<dist{j,1}(PkID{j,1}));
        IDXpksL{j,1}=IDXpksN{j,1}(dist{j,1}(IDXpksN{j,1})>dist{j,1}(PkID{j,1}));            
        
        IDXpksR{j,1}=IDXpksR{j,1}(Int{j,1}(IDXpksR{j,1})>PH(j,1)*PkInt{j,1});
        IDXpksL{j,1}=IDXpksL{j,1}(Int{j,1}(IDXpksL{j,1})>PH(j,1)*PkInt{j,1});
        IDXpksR{j,1}=max(IDXpksR{j,1});
        IDXpksL{j,1}=min(IDXpksL{j,1});
        
        if isempty(IDXpksR{j,1})==1
            IDXpksR{j,1}=1;
        end
        if isempty(IDXpksL{j,1})==1
            IDXpksL{j,1}=length(dist{j,1});
        end
        
        if IDXpksR{j,1}<IDXpksL{j,1}
            if IDXpksR{j,1}>=min(IDXfit{j,1})
                IDXpksR{j,1}=min(IDXval{j,1}(IDXval{j,1}>IDXpksR{j,1}));
                IDXfit{j,1}=[IDXpksR{j,1}:max(IDXfit{j,1})];
            end
            if IDXpksL{j,1}<=max(IDXfit{j,1})
                IDXpksL{j,1}=max(IDXval{j,1}(IDXval{j,1}<IDXpksL{j,1}));
                IDXfit{j,1}=[min(IDXfit{j,1}):IDXpksL{j,1}];
            end
        end
        if isempty(IDXfit{j,1})==1
            if min(SignalIdx{j,1})>1
                minIDXfit=min(SignalIdx{j,1})-1;
            else
                minIDXfit=1;
            end
            if max(SignalIdx{j,1})<length(dist{j,1})
                maxIDXfit=max(SignalIdx{j,1})+1;
            else
                maxIDXfit=length(dist{j,1});
            end
            IDXfit{j,1}=[minIDXfit:maxIDXfit];
        end
        
        % background correction
        IDXfilN{j,1}=intersect(IDXfit{j,1},NoiseIdx{j,1});
        if isempty(NoiseIdx{j,1})==0
            Int{j,1}=Int{j,1}-median(Int{j,1}(NoiseIdx{j,1}));
        else
            Int{j,1}=Int{j,1}-min(Int{j,1}(IDXfit{j,1}));
        end
        
        dist{j,1}=dist{j,1}(IDXfit{j,1});
        Int{j,1}=Int{j,1}(IDXfit{j,1});
        
        % fit with gaussian
        trywhole=1;
        Err=1;
        while Err==1;
            try
                fitG{j,1}=fit(dist{j,1}, (Int{j,1}), 'gauss1');
                Param{j,1}=coeffvalues(fitG{j,1});
                FWHM{j,1}=2.354*(1/sqrt(2))*Param{j,1}(3);
                CentreFit{j,1}=Param{j,1}(2);
%                 Total=Total+1;
                Err=0;

                if FWHM{j,1}>(1.2*(sqrt((ContourSegment(j,1)-ContourSegment(j,3))^2 + (ContourSegment(j,2)-ContourSegment(j,4))^2)))
                    ME=error('too big');
                    Err=1;
                end
                    
            catch ME
%                 display(ME);
                    if trywhole==1
                            Err=1;
                            trywhole=0;
                            dist{j,1}=distO{j,1};
                            Int{j,1}=c{j,1};
                            if isempty(NoiseIdx{j,1})==0
                                Int{j,1}=Int{j,1}-median(Int{j,1}(NoiseIdx{j,1}));
                            else
                                Int{j,1}=Int{j,1}-min(Int{j,1});
                            end
                    else
                        messedupframes=[messedupframes j];
                         Err=0;
%                          display(sprintf('Error in frame %d segment %d: please select profile', i,j));
%                          Fig=figure;
%                          subplot(1,2,1)
%                          imshow(OriginalImage,[])
%                          hold on
%                          line([ContourSegment(j,1), ContourSegment(j,3)], [ContourSegment(j,2),ContourSegment(j,4)]);
%                          subplot(1,2,2)
%                          plot(distO{j,1},(c{j,1}-min(c{j,1})))
%                          title('Select left cutoff position')
%                         %      datacursormode on
%                         %      dcm_obj = datacursormode;
%                         [xL,yL]=ginput(1);
%                         title('Select right cutoff position')
%                         [xR yR]=ginput(1);
%                         IntO{j,1}=(c{j,1}-min(c{j,1}));
%                         dist{j,1}=distO{j,1};
%                         IDX2{j,1}=[1:length(dist{j,1})];
%                         IDX2{j,1}(dist{j,1}<xL)=[];
%                         Int{j,1}(dist{j,1}<xL)=[];
%                         dist{j,1}(dist{j,1}<xL)=[];
%                         Int{j,1}(dist{j,1}>xR)=[];
%                         IDX2{j,1}(dist{j,1}>xR)=[];
%                         dist{j,1}(dist{j,1}>xR)=[];
%                         
%                         if isempty(NoiseIdx{j,1})==0
%                             Int{j,1}=Int{j,1}-median(IntO{j,1}(NoiseIdx{j,1}));
%                         else
%                             Int{j,1}=Int{j,1}-min(Int{j,1});
%                         end
%                         close(Fig)
                        clear ME
                    end
            end
        end
    end
% end
% close(h)

end