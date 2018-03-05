function [TubeRnew, EnvRnew, EnvBackbone] = MinECurvatureExtract(segCurv1, segCurv2, StartFrame, EndFrame, segLength, smallestID, segmentRadii, Extents, k1bb)
%outputs: min envelope curvature at min diameter (MinECconSite) 
%minnenvelope curvature at min diameter averaged with 2 nearest points (smoothMinECconSite) 
% average of envelope curvatures in segment corresponding to min D (  avMinECconSite)

% Totlength=segLength{1,1}(smallestID)+segLength{2,1}(smallestID);
% WeightS1=segLength{1,1}(smallestID)/Totlength;
% WeightS2=segLength{2,1}(smallestID)/Totlength;

for f=StartFrame:EndFrame
    
%         TubeR(f, :)=segmentRadii{1,f}(smallestID(f));
        
    
    
    for len=1:length(Extents{1,f})
    
        Totlength(f,len)=sum(segLength{1,f}(Extents{1,f}(len, 1):Extents{1,f}(len, 2)))+sum(segLength{2,f}(Extents{1,f}(len, 1):Extents{1,f}(len, 2)));
        WeightS1(f,len)=sum(segLength{1,f}(Extents{1,f}(len, 1):Extents{1,f}(len, 2)))/Totlength(f,len);
        WeightS2(f,len)=sum(segLength{2,f}(Extents{1,f}(len, 1):Extents{1,f}(len, 2)))/Totlength(f,len);
        TubeRnew(f, len)=mean(segmentRadii{1,f}(Extents{1,f}(len, 1):Extents{1,f}(len, 2))); 

%         MinECconSite(f,1)= min(segCurv1{f,1});
%         MinECconSite(f,2)= min(segCurv2{f,1});
        
    avMinECconSite1(f,len)=mean(segCurv1{f,len});  %sum(segCurv1{f,len})/size(segCurv1{f,len},2);
    avMinECconSite2(f,len)=mean(segCurv2{f,len});  %sum(segCurv2{f,len})/size(segCurv2{f,len},2);
 
        avMinECconSiteTot(f,len)= WeightS1(f,len)*avMinECconSite1(f,len)+ WeightS2(f,len)*avMinECconSite2(f,len);
        EnvRnew(f,len)=1./avMinECconSiteTot(f,len);
        
        
        %%envelope curvature of backbone
        EnvBackbone(f,len)=mean(k1bb{1,f}(Extents{1,f}(len, 1):Extents{1,f}(len, 2)));
        
end
end

% 
% %take average of the min curvature point with its 2 nearest neighbours (1 to the left and the other to the right)
% for f=StartFrame:EndFrame
%         minTag(f,1)=find(segCurv1{f,1}==MinECconSite(f,1));
%         minTag(f,2)=find(segCurv2{f,1}==MinECconSite(f,2));
% end
% 
% smoothMinECconSite=[];
% for f=StartFrame:EndFrame
%        %side1
%         %if minimum curvature found at edge, take its 2 nearest points to
%         %the left or right
%         if minTag(f,1)==1&& size(segCurv1{f,1},2)>=3
%         smoothMinECconSite(f,1)=sum(segCurv1{f,1}(1,1:3))/3;
%         end
%         
%        if minTag(f,1)==size(segCurv1{f,1},2)&& size(segCurv1{f,1},2)>=3
%         smoothMinECconSite(f,1)=sum(segCurv1{f,1}(1,end-2:end))/3;
%        end
%         
%         %otherwise, average the min with a point to its left and another to
%         %its right
%           if minTag(f,1)~=1 && minTag(f,1)~=size(segCurv1{f,1},2)&& size(segCurv1{f,1},2)>=3
%            smoothMinECconSite(f,1)=sum(segCurv1{f,1}(1,minTag(f,1)-1:minTag(f,1)+1))/3;
%           end
%         
%           if size(segCurv1{f,1},2)<3
%               smoothMinECconSite(f,1)=MinECconSite(f,1);
%           end
%         
%         %side 2
%          %if minimum curvature found at edge, take its 2 nearest points to
%         %the left or right
%         if minTag(f,2)==1&& size(segCurv2{f,1},2)>=3
%         smoothMinECconSite(f,2)=sum(segCurv2{f,1}(1,1:3))/3;
%         end
%         if minTag(f,2)==size(segCurv2{f,1},2)&& size(segCurv2{f,1},2)>=3
%         smoothMinECconSite(f,2)=sum(segCurv2{f,1}(1,end-2:end))/3;
%         end
%         
%         %otherwise, average the min with a point to its left and another to
%         %its right
%         
%         if minTag(f,2)~=1 && minTag(f,2)~=size(segCurv2{f,1},2)&& size(segCurv2{f,1},2)>=3
%            smoothMinECconSite(f,2)=sum(segCurv2{f,1}(1,minTag(f,2)-1:minTag(f,2)+1))/3;
%         end
%         
%         if size(segCurv2{f,1},2)<3
%               smoothMinECconSite(f,2)=MinECconSite(f,2);
%           end
%               
% end
% 
% 
% 
