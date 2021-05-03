%% 实验函数
 % 点云数据为 3 组 
 
%  归一化配准实验
%  20210501
function GC_test( vase_data1 , DTran_vase_data2 , DTran_vase_data3 , correspondence21 , correspondence31 , correspondence32)

%   data1 L1 -> data2 L3
%   data1 L1 -> data3 L4
%   data2 L1 -> data3 L3

%% 计算初始平移系数

T2 = eye(4);
T3 = eye(4);

reg12_data1_L1 = vase_data1.SubBoundaryInformation.Points(1).Line;
reg12_data2_L3 = DTran_vase_data2.SubBoundaryInformation.Points(3).Line;
len_reg12_data2_L3 = length(correspondence21);

reg13_data1_L1 = vase_data1.SubBoundaryInformation.Points(1).Line;
reg13_data3_L4 = DTran_vase_data3.SubBoundaryInformation.Points(4).Line;
len_reg13_data3_L4 = length(correspondence31);

reg23_data2_L1 = DTran_vase_data2.SubBoundaryInformation.Points(1).Line;
reg23_data3_L3 = DTran_vase_data3.SubBoundaryInformation.Points(3).Line;
len_reg23_data3_L3 = length(correspondence32);


totolnum = len_reg12_data2_L3 + len_reg13_data3_L4 + len_reg23_data3_L3 ;

%% KDTREE

reg13_data1_L1_kdtreeobj = KDTreeSearcher( reg13_data1_L1 , 'distance' , 'euclidean' );    %建立 Refercloud的kdtree  
reg23_data2_L1_kdtreeobj = KDTreeSearcher( reg23_data2_L1 , 'distance' , 'euclidean' );    %建立 Refercloud的kdtree

%% SHOW
% pcshow(reg23_data3_L3,[1,0,0],'MarkerSize',100)
% hold on
% % pcshow( vase_data1.pc.Location , [0,0,0] , 'MarkerSize' , 10 );
% % pcshow( DTran_vase_data2.Transed , [0,0,0] , 'MarkerSize' , 10 );
% % pcshow( DTran_vase_data3.Transed , [0,0,0] , 'MarkerSize' , 10 );
% pcshow(reg23_data2_L1,[1,0,0],'MarkerSize',100)
% pcshow(reg12_data1_L1(correspondence21,1:3),[0,1,0],'MarkerSize',100)
% pcshow(reg13_data1_L1(correspondence31,1:3),[0,0,1],'MarkerSize',100)
% pcshow(reg12_data2_L3,[0,1,0],'MarkerSize',100)
% pcshow(reg13_data3_L4,[0,0,1],'MarkerSize',100)
% hold off

% max = 100;

       
    %% 迭代计算
    
    Iterativetimes_num = 1;
    %     H = zeros(6*3);  % 定义 Hessian 矩阵
    b = ones( ( len_reg12_data2_L3 + len_reg13_data3_L4 + len_reg23_data3_L3 ) * 3 , 1 ); % 定义 J * e
    g = b;
    while( norm(g) > 0.000006 )   %迭代计算开始
        
%         if(Iterativetimes_num > 500)
%             
%             return;
%             
%         end
                
%     H = zeros( 6*3 );  % 定义 Hessian 矩阵
    b = zeros( ( len_reg12_data2_L3 + len_reg13_data3_L4 + len_reg23_data3_L3 ) * 3 , 1); % 定义 J * e 
    J = zeros( ( len_reg12_data2_L3 + len_reg13_data3_L4 + len_reg23_data3_L3 ) * 3 , 2 * 6 );
    
    %% 构建 J b
    
    for i = 1 : len_reg12_data2_L3
        
        J( ((i-1)*3 + 1):(i*3) , 1:6 ) = [ eye(3) , - ToAntisymmetric_Mat( reg12_data2_L3( i , 1:3 )' ) ];
        b( ((i-1)*3 + 1):(i*3) , 1 ) = ( reg12_data2_L3(i,1:3)' - reg12_data1_L1( correspondence21(i,1) , 1:3)' );
        
    end
    
    j = 1;
    for i = (len_reg12_data2_L3 + 1) : (len_reg12_data2_L3 + len_reg13_data3_L4)
        
        J( ((i-1)*3 + 1):(i*3) , 7:12 ) = [ eye(3) , - ToAntisymmetric_Mat( reg13_data3_L4( j , 1:3 )' ) ];
        b( ((i-1)*3 + 1):(i*3) , 1 ) = ( reg13_data3_L4(j,1:3)' - reg13_data1_L1( correspondence31(j,1) , 1:3)' );
        j = j + 1;
        
    end
    
    j = 1;
    for i = (len_reg12_data2_L3 + len_reg13_data3_L4 + 1) : (totolnum)
        
        J( ((i-1)*3 + 1):(i*3) , 1:6 ) = [ -eye(3) ,  ToAntisymmetric_Mat( reg23_data2_L1( j , 1:3 )' ) ];
        J( ((i-1)*3 + 1):(i*3) , 7:12 ) = [ eye(3) , - ToAntisymmetric_Mat( reg23_data3_L3( j , 1:3 )' ) ];
        
        b( ((i-1)*3 + 1):(i*3) , 1 ) = ( reg23_data3_L3(j,1:3)' - reg23_data2_L1( correspondence32(j,1) , 1:3)' );
        j = j + 1;
        
    end
    
  
    
%% H x = -g
     
        
    H = J' * J;   % H 海瑟矩阵
    g = J' * b;   % g = J' * b
    
    
    delt = H \ -g;
    
    T2_so = delt( 1 : 6  , 1 );
    T3_so = delt( 7 : 12 , 1 );
    
    T2_deltT = so2SO(T2_so);
    T3_deltT = so2SO(T3_so);
    
    T2 = T2_deltT * T2;
    T3 = T3_deltT * T3;
       
    reg12_data2_L3 = Trans( reg12_data2_L3 , T2_deltT );
    reg13_data3_L4 = Trans( reg13_data3_L4 , T3_deltT );
    reg23_data2_L1 = Trans( reg23_data2_L1 , T2_deltT );
    reg23_data3_L3 = Trans( reg23_data3_L3 , T3_deltT );
    
    pcshow(reg23_data3_L3,[1,0,0],'MarkerSize',100)
    hold on
    
%     pcshow( vase_data1.pc.Location , [0,0,0] , 'MarkerSize' , 10 );
%     pcshow( DTran_vase_data2.Transed , [0,0,0] , 'MarkerSize' , 10 );
%     pcshow( DTran_vase_data3.Transed , [0,0,0] , 'MarkerSize' , 10 );
    
    pcshow(reg23_data2_L1,[1,0,0],'MarkerSize',100)
    pcshow(reg12_data1_L1(correspondence21,1:3),[0,1,0],'MarkerSize',100)
    pcshow(reg13_data1_L1(correspondence31,1:3),[0,0,1],'MarkerSize',100)
    pcshow(reg12_data2_L3,[0,1,0],'MarkerSize',100)
    pcshow(reg13_data3_L4,[0,0,1],'MarkerSize',100)
    hold off
    
    Iterativetimes_num = Iterativetimes_num + 1;

%     if norm(g) < max
%         
%         max = norm(g)
%         
%     end
     norm(g)
     norm(b)
     
     for i = 1 : length(reg13_data3_L4)   %确定点对应关系
         
         [ idx , ~ ] = knnsearch( reg13_data1_L1_kdtreeobj , reg13_data3_L4(i,1:3) , 'dist' , 'euclidean' , 'k' , 1);
         correspondence31(i,1) = idx;
         
     end
     
%      for i = 1 : length(reg23_data3_L3)   %确定点对应关系
%          
%          [ idx , ~ ] = knnsearch( reg23_data2_L1_kdtreeobj , reg23_data3_L3(i,1:3) , 'dist' , 'euclidean' , 'k' , 1);
%          correspondence32(i,1) = idx;
%          
%      end
     
     
     
     
    end
    
    
    



end


function [deltT] = so2SO(so)

% delt = H \ -b;

t_ = so(1:3,1);

fa = so(4:6,1);

norm_delt = norm(fa);
a = fa / norm_delt;

J_ = ( sin( norm_delt ) / norm_delt ) * eye(3) + ( 1 - sin(norm_delt) / norm_delt )*( a * a' )+(( 1 - cos(norm_delt) ) / norm_delt ) * ToAntisymmetric_Mat(a);
deltR = cos( norm_delt ) * eye(3) + ( 1 - cos( norm_delt ) ) * ( a * a' ) + sin( norm_delt ) * ToAntisymmetric_Mat(a);


deltT = [deltR     ,  J_ * t_;
    zeros(1,3) , 1       ];


end



function [cloud] = Trans(cloud,T)

    [l , w] = size(cloud);
    
    if ( w ~= 3 )
        
        cloud =cloud';
        [l , ~] = size(cloud);
        
    end
    
    
    cloud(:,4) = ones(l,1);
    
    cloud = ( T * cloud' )';
    
    cloud = cloud( :, 1:3 );



end




function [Mat] = ToAntisymmetric_Mat(vector)
    
    [ w , ~ ] = size(vector);
    
    if(w ~= 3)
        
        vector = vector';
        
    end
       
    a1 = vector(1,1);
    a2 = vector(2,1);
    a3 = vector(3,1);
    
    
    Mat = [ 0  , -a3 ,  a2;
            a3 ,  0  , -a1;
           -a2 ,  a1 ,  0];

end



 