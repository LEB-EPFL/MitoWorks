function BreakingFrame=splitLastFrame(SegmentedMitoIDImage,BreakFrame,IDmin,subNormals)

    imshow((SegmentedMitoIDImage{1,BreakFrame}));
    line(subNormals{1,BreakFrame}(IDmin(BreakFrame,1),1:2:3),subNormals{1,BreakFrame}(IDmin(BreakFrame,1),2:2:4), 'LineWidth', 2, 'Color', 'k');
    
    BreakingFrame=im2bw(frame2im((getframe)),graythresh(frame2im((getframe))));
    
    BreakingFrame(:,1:2)=0;
    BreakingFrame(:,(end-1):end)=0;
    BreakingFrame(1:2,:) = 0;
    BreakingFrame((end-1):end,:)=0;
%     size(SegmentedMitoIDImage{1,BreakFrame},1:2)
x=size(SegmentedMitoIDImage{1,BreakFrame});
    BreakingFrame = imresize(logical(BreakingFrame),x(1:2));
    
end
