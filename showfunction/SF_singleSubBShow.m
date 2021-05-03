function SF_singleSubBShow(pc,SubB_indx,Normaltag)
    

    pcshow(pc.pc.Location,[0,0,0],'MarkerSize',10 );
    hold on
    pcshow( pc.SubBoundaryInformation.Points(SubB_indx).Line , [1,0,0] , 'MarkerSize' , 100 );
    
    if( Normaltag == true )
        
        Line1 = pc.SubBoundaryInformation.Points(SubB_indx).Line;
        normal1 = pc.SubBoundaryInformation.Points(SubB_indx).Normal;
        
        quiver3(Line1(:,1),Line1(:,2),Line1(:,3),normal1 (:,1),normal1 (:,2),normal1(:,3));
        
    end
    
    hold off


end

