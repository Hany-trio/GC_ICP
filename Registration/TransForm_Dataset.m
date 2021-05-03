%%  配准数据
 % input：objectdata（数据集）  R（translate parameter）
 
 function Tran_pc = TransForm_Dataset( pc , T )
    
    Tran_pc = pc;
    Tran_pc.Transed = Transform_cp( pc.pc.Location , T );
    Tran_pc.Boundaryinformation.BoundaryPoints = Transform_cp( pc.Boundaryinformation.BoundaryPoints , T );
    Tran_pc.TurningPoint.Points = Transform_cp( pc.TurningPoint.Points , T );
    
    [ ~ , l ] = size(pc.SubBoundaryInformation.Points);
    
    
    for i = 1 : l
       
        Tran_pc.SubBoundaryInformation.Points(i).TrPoint = Transform_cp( pc.SubBoundaryInformation.Points(i).TrPoint,T);
        Tran_pc.SubBoundaryInformation.Points(i).Line  = Transform_cp( pc.SubBoundaryInformation.Points(i).Line , T );
        Tran_pc.SubBoundaryInformation.Points(i).Normal  = Transform_cp( pc.SubBoundaryInformation.Points(i).Normal , T );
        
    end
    


end



function p = Transform_cp(p,T)     % 旋转一个维度为3的点集/法线集
       
	[l,w] = size(p);
        
	if  ( l ~= 3 )
           
        p  = p';
        [ ~ , w ] = size( p );
            
	end
        
	p = [ p ; ones(1,w) ];
        
	p = ( T * p )';
        
	p = p(:,1:3);

    
end