clc
clear

%% 加载数据 

  % 该实验主要是三组点云数据
  % 加载的三组数据非 .pcd 文件，而是特征提取后的数据
  
    vase_data1 = load('C:\Users\HanY、\Desktop\Global consistent\data\featuredata\vase_sample2-Dousubsample20_1');
    Tran_vase_data2 = load('C:\Users\HanY、\Desktop\Global consistent\data\featuredata\tran_vase_sample2-Dousubsample20_2');
    Tran_vase_data3 = load('C:\Users\HanY、\Desktop\Global consistent\data\featuredata\tran_vase_sample2-Dousubsample20_3');
  
    vase_data1 = vase_data1.pc_FeatureInformation_;
    Tran_vase_data2 = Tran_vase_data2.pc_FeatureInformation_;
    Tran_vase_data3 = Tran_vase_data3.pc_FeatureInformation_;
  
%% 计算初始位姿参数
    
    %可视化函数  
    %   data1 L1 -> data2 L3
    %   data1 L1 -> data3 L4
    %   data2 L1 -> data3 L3
    
    
%     SF_singleSubBShow(vase_data1,1,false);
%     hold on
%     SF_singleSubBShow(Tran_vase_data2,1,false);
%     hold on
%     SF_singleSubBShow(Tran_vase_data3,3,false);
    
    % [Tij] represent the tranformation parameter of object i to object j; i -> j;  
    % [correspondenceij] represent the correspondence points of data j to
    %   data i, length(correspondenceij) = length(dataj_registLine_)
    [ correspondence21 , T21 ] = GC_BaseLine_Register( vase_data1 , Tran_vase_data2 , 1 , 3 , "NICP" );   % 初始参数计算 
    [ correspondence31 , T31 ] = GC_BaseLine_Register( vase_data1 , Tran_vase_data3 , 1 , 4 , "NICP" );
    [ correspondence32 , T32 ] = GC_BaseLine_Register( Tran_vase_data2 , Tran_vase_data3 , 1 , 3 , "NICP" );
    %     T21 和 T31 分别是Tran_vase_data2和Tran_vase_data3 的 初始变换参数
    
    
%% 全局一致性估计
    
DTran_vase_data2 = TransForm_Dataset( Tran_vase_data2 , T21 );
DTran_vase_data3 = TransForm_Dataset( Tran_vase_data3 , T31 );

% 1
% 
% pcshow( vase_data1.pc.Location , [1,0,0] , 'MarkerSize' , 50 );
% hold on
% pcshow( DTran_vase_data2.Transed , [0,1,0] , 'MarkerSize' , 50 );
% pcshow( DTran_vase_data3.Transed , [0,0,1] , 'MarkerSize' , 50 );


% 

GC_test( vase_data1 , DTran_vase_data2 , DTran_vase_data3 , correspondence21 , correspondence31 , correspondence32)

% 已知 变换 初始化参数 T21 T31
% 已知 对应点 初始化参数 correspondence21 correspondence31 correspondence32
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  