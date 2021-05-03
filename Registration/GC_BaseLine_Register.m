function [ correspondence , T ] = GC_BaseLine_Register( pc1 , pc2 , pc1_Lindx , pc2_Lindx , Registag )
    
    Refer_Line = pc1.SubBoundaryInformation.Points(pc1_Lindx).Line;
    Match_Line = pc2.SubBoundaryInformation.Points(pc2_Lindx).Line;

    Refer_Line_Normal = pc1.SubBoundaryInformation.Points(pc1_Lindx).Normal;
    Match_Line_Normal = pc2.SubBoundaryInformation.Points(pc2_Lindx).Normal;

    if(Registag == "ICP")
       
        [ correspondence , T ] = n_ICP( Refer_Line , Match_Line );
        
    elseif( Registag == "NICP" )
            
        [ correspondence , T ] = HCode_NICP( Refer_Line , Match_Line , Refer_Line_Normal , Match_Line_Normal , 0.003 );
        
    else
        
        [T] = eye(4);
        
    end
    
    
    
end

