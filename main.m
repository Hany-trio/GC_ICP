clc
clear

%% �������� 

  % ��ʵ����Ҫ�������������
  % ���ص��������ݷ� .pcd �ļ�������������ȡ�������
  
    vase_data1 = load('C:\Users\HanY��\Desktop\Global consistent\data\featuredata\vase_sample2-Dousubsample20_1');
    Tran_vase_data2 = load('C:\Users\HanY��\Desktop\Global consistent\data\featuredata\tran_vase_sample2-Dousubsample20_2');
    Tran_vase_data3 = load('C:\Users\HanY��\Desktop\Global consistent\data\featuredata\tran_vase_sample2-Dousubsample20_3');
  
    vase_data1 = vase_data1.pc_FeatureInformation_;
    Tran_vase_data2 = Tran_vase_data2.pc_FeatureInformation_;
    Tran_vase_data3 = Tran_vase_data3.pc_FeatureInformation_;
  
%% �����ʼλ�˲���
    
    %���ӻ�����  
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
    [ correspondence21 , T21 ] = GC_BaseLine_Register( vase_data1 , Tran_vase_data2 , 1 , 3 , "NICP" );   % ��ʼ�������� 
    [ correspondence31 , T31 ] = GC_BaseLine_Register( vase_data1 , Tran_vase_data3 , 1 , 4 , "NICP" );
    [ correspondence32 , T32 ] = GC_BaseLine_Register( Tran_vase_data2 , Tran_vase_data3 , 1 , 3 , "NICP" );
    %     T21 �� T31 �ֱ���Tran_vase_data2��Tran_vase_data3 �� ��ʼ�任����
    
    
%% ȫ��һ���Թ���
    
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

% ��֪ �任 ��ʼ������ T21 T31
% ��֪ ��Ӧ�� ��ʼ������ correspondence21 correspondence31 correspondence32
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  