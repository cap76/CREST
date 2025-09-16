%addpath(genpath('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/code/gpml-matlab-v3.6-2015-07-07'))
cd('/Users/christopherpenfold/Desktop/Code/Spatial modelling/SpatialModelling/Functions')
cd('/Users/christopherpenfold/Desktop/Code/Spatial modelling/SpatialModelling')
addpath(genpath('/Users/christopherpenfold/Desktop/Code/Spatial\ modelling/SpatialModelling/'))

cd('/Users/christopherpenfold/Desktop/Code/Spatial modelling/SpatialModelling/Functions')
addpath(genpath('./'))
addpath(genpath('/Users/christopherpenfold/Desktop/Code/Spatial modelling/SpatialModelling'))
cd('/Users/christopherpenfold/Desktop/Code/Spatial modelling/SpatialModelling')



%Cols =  c("EmDisc"="#2D7EE0",
%          "Am/EmDisc"="#00B9E3",
%          "Am"="#3CDDDD",
%          "PGC"="#E154E8",
%          "Hypoblast"="#FF8C64",
%          "YS"="#2E892E",
%          "ExMes/YS"="#22BF35",
%          "ExMes"="#9FD13F",
%          "CTB"="#E2C120",
%          "STB"="#D68C34",
%          "EVT"="#BF3737",
%          "Stromal"="#AD71BC",
%          "Epithelial"="#C4807F",
%          "Ciliated"="#7F3C3C")

%Load CS6
[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,ShotsCS6] = LoadShots('CS6');

File1 = '/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping5/S1_set1_all_withPreImpl_byCl.csv'
File2 = '/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping5/C1_set1_withPreImpl_byCl.csv'
[newX1,newY1,newZ1,X1,Y1,Z1,W1,CellType1,CellTypeReg] = ProjectData(File1,File2,ShotsCS6,10);

%csvwrite('~/Desktop/CT.csv',CellType1)
%csvwrite('~/Desktop/CTR.csv',CellTypeReg)
%writecell(CellType1,'~/Desktop/CT.csv') 
%writecell(CellTypeReg,'~/Desktop/CTR.csv') 

l0 = find( isnan(newX1)==0 );


S1 = readtable(File1);
C1 = readtable(File2);


%l1 = find( isnan(newX1)==0 & (strcmp(CellType1,'Am')==1 ) );
%l2 = find( isnan(newX1)==0 & (strcmp(CellType1,'Am/EmDisc')==1 ) );
%l3 = find( isnan(newX1)==0 & (strcmp(CellType1,'EmDisc')==1 ) );
%l4 = find( isnan(newX1)==0 & (strcmp(CellType1,'ExMes')==1 ) );
%l5 = find( isnan(newX1)==0 & (strcmp(CellType1,'ExMes/SYS')==1 ) );
%l6 = find( isnan(newX1)==0 & (strcmp(CellType1,'Hypoblast')==1 ) );
%l7 = find( isnan(newX1)==0 & (strcmp(CellType1,'PGC')==1 ) );
%l8 = find( isnan(newX1)==0 & (strcmp(CellType1,'SYS')==1 ) );


l1 = find( isnan(newX1)==0 & (strcmp(CellType1,'Amnion')==1 ) );
l2 = find( isnan(newX1)==0 & (strcmp(CellType1,'EmbDisc')==1 ) );
l3 = find( isnan(newX1)==0 & (strcmp(CellType1,'ExEmbMes')==1 ) );
l4 = find( isnan(newX1)==0 & (strcmp(CellType1,'Hypoblast')==1 ) );
l5 = find( isnan(newX1)==0 & (strcmp(CellType1,'YolkSac')==1 ) );

l6 = find( isnan(newX1)==0 & (strcmp(CellType1,'STB')==1 ) );
l7 = find( isnan(newX1)==0 & (strcmp(CellType1,'CTB')==1 ) );
l8 = find( isnan(newX1)==0 & (strcmp(CellType1,'eCTB')==1 ) );



%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[113, 218, 219]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[45, 122, 234]/255,'LineStyle','none','FaceAlpha',1); 
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[169, 208, 88]/255,'LineStyle','none','FaceAlpha',1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[240, 146, 109]/255,'LineStyle','none','FaceAlpha',1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[73, 135, 59]/255,'LineStyle','none','FaceAlpha',1);
f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[203, 143, 70]/255,'LineStyle','none','FaceAlpha',1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
colormap(parula);
print(h,['~/Desktop/Am_CS6_AllTissues_withXiang_background.png'],'-dpng','-r600')



%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[113, 218, 219]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[45, 122, 234]/255,'LineStyle','none','FaceAlpha',1); 
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[169, 208, 88]/255,'LineStyle','none','FaceAlpha',1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[240, 146, 109]/255,'LineStyle','none','FaceAlpha',1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[73, 135, 59]/255,'LineStyle','none','FaceAlpha',1);
f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[203, 143, 70]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
colormap(parula);
print(h,['~/Desktop/Am_CS6_AllTissues_withXiang_backgroundB.png'],'-dpng','-r600')


%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[113, 218, 219]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[45, 122, 234]/255,'LineStyle','none','FaceAlpha',.1); 
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[169, 208, 88]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[240, 146, 109]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[73, 135, 59]/255,'LineStyle','none','FaceAlpha',.1);
f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[203, 143, 70]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
%scatter3(newX1(l0,1) + rand(length(l0),1)*1,newY1(l0,1)+ rand(length(l0),1)*1,newZ1(l0,1)+ rand(length(l0),1)*1,W1(l0)*60,[12,156,245]/255,'filled')

scatter3(newX1(l1,1) + rand(length(l1),1)*1,newY1(l1,1)+ rand(length(l1),1)*1,newZ1(l1,1)+ rand(length(l1),1)*1,W1(l1)*60,[60, 221, 221]/255,'filled')
scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[0, 185, 227]/255,'filled')
scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[159, 209, 63]/255,'filled')
scatter3(newX1(l4,1) + rand(length(l4),1)*1,newY1(l4,1)+ rand(length(l4),1)*1,newZ1(l4,1)+ rand(length(l4),1)*1,W1(l4)*60,[255, 140, 100]/255,'filled')
scatter3(newX1(l5,1) + rand(length(l5),1)*1,newY1(l5,1)+ rand(length(l5),1)*1,newZ1(l5,1)+ rand(length(l5),1)*1,W1(l5)*60,[46, 137, 46]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[117,76,36]/255,'filled')

scatter3(newX1(l6,1) + rand(length(l6),1)*1,newY1(l6,1)+ rand(length(l6),1)*1,newZ1(l6,1)+ rand(length(l6),1)*1,W1(l6)*60,[214, 140, 52]/255,'filled')
scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,[226, 193, 32]/255,'filled')
scatter3(newX1(l8,1) + rand(length(l8),1)*1,newY1(l8,1)+ rand(length(l8),1)*1,newZ1(l8,1)+ rand(length(l8),1)*1,W1(l8)*60,[226, 193, 32]/255,'filled')


colormap(parula);
print(h,['~/Desktop/Am_CS6_AllTissues_withXiang2.png'],'-dpng','-r600')

print(h,['~/Desktop/Am_CS6_AllTissues_withXiang3.png'],'-dpng','-r600')

% 
% D = zeros(length(l0),8);
% 
% for i = 1:length(l0)
%     llocal = l0(i);
% 
%    x = OBJ2.vertices(OBJ2.objects(20).data.vertices,1);
%    dis(i,1)= min( (x - repmat(newX1(llocal,1),size(x,1),1)).^2 + ( OBJ2.vertices(OBJ2.objects(20).data.vertices,2) - repmat(newY1(llocal,1),size(x,1),1)).^2  + ( OBJ2.vertices(OBJ2.objects(20).data.vertices,3) - repmat(newZ1(llocal,1),size(x,1),1)).^2  );
% 
%    x = OBJ2.vertices(OBJ2.objects(4).data.vertices,1);
%    dis(i,2)= min( (x - repmat(newX1(llocal,1),size(x,1),1)).^2 + ( OBJ2.vertices(OBJ2.objects(4).data.vertices,2) - repmat(newY1(llocal,1),size(x,1),1)).^2  + ( OBJ2.vertices(OBJ2.objects(4).data.vertices,3) - repmat(newZ1(llocal,1),size(x,1),1)).^2  );
% 
%    x = OBJ2.vertices(OBJ2.objects(8).data.vertices,1);
%    dis(i,3)= min( (x - repmat(newX1(llocal,1),size(x,1),1)).^2 + ( OBJ2.vertices(OBJ2.objects(8).data.vertices,2) - repmat(newY1(llocal,1),size(x,1),1)).^2  + ( OBJ2.vertices(OBJ2.objects(8).data.vertices,3) - repmat(newZ1(llocal,1),size(x,1),1)).^2  );
% 
%    x = OBJ2.vertices(OBJ2.objects(24).data.vertices,1);
%    dis(i,4)= min( (x - repmat(newX1(llocal,1),size(x,1),1)).^2 + ( OBJ2.vertices(OBJ2.objects(24).data.vertices,2) - repmat(newY1(llocal,1),size(x,1),1)).^2  + ( OBJ2.vertices(OBJ2.objects(24).data.vertices,3) - repmat(newZ1(llocal,1),size(x,1),1)).^2  );
% 
%    x = OBJ2.vertices(OBJ2.objects(16).data.vertices,1);
%    dis(i,5)= min( (x - repmat(newX1(llocal,1),size(x,1),1)).^2 + ( OBJ2.vertices(OBJ2.objects(16).data.vertices,2) - repmat(newY1(llocal,1),size(x,1),1)).^2  + ( OBJ2.vertices(OBJ2.objects(16).data.vertices,3) - repmat(newZ1(llocal,1),size(x,1),1)).^2  );
% 
%    x = OBJ2.vertices(OBJ2.objects(12).data.vertices,1);
%    dis(i,6)= min( (x - repmat(newX1(llocal,1),size(x,1),1)).^2 + ( OBJ2.vertices(OBJ2.objects(12).data.vertices,2) - repmat(newY1(llocal,1),size(x,1),1)).^2  + ( OBJ2.vertices(OBJ2.objects(12).data.vertices,3) - repmat(newZ1(llocal,1),size(x,1),1)).^2  );
% 
%    x = OBJ2.vertices(OBJ2.objects(28).data.vertices,1);
%    dis(i,7)= min( (x - repmat(newX1(llocal,1),size(x,1),1)).^2 + ( OBJ2.vertices(OBJ2.objects(28).data.vertices,2) - repmat(newY1(llocal,1),size(x,1),1)).^2  + ( OBJ2.vertices(OBJ2.objects(28).data.vertices,3) - repmat(newZ1(llocal,1),size(x,1),1)).^2  );
% 
%   x = OBJ2.vertices(OBJ2.objects(32).data.vertices,1);
%    dis(i,8)= min( (x - repmat(newX1(llocal,1),size(x,1),1)).^2 + ( OBJ2.vertices(OBJ2.objects(32).data.vertices,2) - repmat(newY1(llocal,1),size(x,1),1)).^2  + ( OBJ2.vertices(OBJ2.objects(32).data.vertices,3) - repmat(newZ1(llocal,1),size(x,1),1)).^2  );
% 
% 
% end
% 
% 
% writecell(CellType1,'~/Desktop/CT_CS6.csv') 
% writecell(CellTypeReg,'~/Desktop/CTR_CS6.csv') 
% writematrix(dis,'~/Desktop/Distances_CS6.csv') 
% writetable(C1(:,1),'~/Desktop/Cell_CS6.csv') 
% writematrix(dis,'~/Desktop/Distances_CS6.csv') 
% 



%And direct mapping. This looks horrible in the section. Best to do this as
%a 3D thing, I think
[OBJ1,section] = LoadCS5('3D');
[D,Locations,XYZ,CellType,ShotsCS5] = LoadShots('CS5');
[OutputCS5] = loadCS5Scaffold(D,Locations,ShotsCS5);

File1 = '/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping5/S1_set1_all_withPreImpl_byCl.csv'
File2 = '/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping5/C1_set1_withPreImpl_byCl.csv'
[newX1,newY1,newZ1,X1,Y1,Z1,W1,CellType1] = ProjectDataCS5(File1,File2,ShotsCS5,10);



l1 = find( isnan(newX1)==0 & (strcmp(CellType1,'Amnion')==1 ) );
l2 = find( isnan(newX1)==0 & (strcmp(CellType1,'EmbDisc')==1 ) );
l3 = find( isnan(newX1)==0 & (strcmp(CellType1,'ExEmbMes')==1 ) );
l4 = find( isnan(newX1)==0 & (strcmp(CellType1,'Hypoblast')==1 ) );
l5 = find( isnan(newX1)==0 & (strcmp(CellType1,'YolkSac')==1 ) );
l6 = find( isnan(newX1)==0 & (strcmp(CellType1,'STB')==1 ) );
l7 = find( isnan(newX1)==0 & (strcmp(CellType1,'CTB')==1 ) );
l8 = find( isnan(newX1)==0 & (strcmp(CellType1,'eCTB')==1 ) );

opac = 0.1
l0 = find( isnan(newX1)==0 );



h = figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[113, 218, 219]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[0, 185, 227]/255,'LineStyle','none','FaceAlpha',1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240, 146, 109]/255,'LineStyle','none','FaceAlpha',1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[73, 135, 59]/255,'LineStyle','none','FaceAlpha',1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[169, 208, 88]/255,'LineStyle','none','FaceAlpha',1);
f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[203, 143, 70]/255,'LineStyle','none','FaceAlpha',1);

axis equal
axis off
%view([a1 b1])
view([[11.1385,-68.9693]])
camlight('left')
material dull 

colormap(parula);
print(h,['~/Desktop/Am_CS5_AllTissues_withXiangBckGr.png'],'-dpng','-r600')





h = figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[113, 218, 219]/255,'LineStyle','none','FaceAlpha',opac);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[0, 185, 227]/255,'LineStyle','none','FaceAlpha',opac);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240, 146, 109]/255,'LineStyle','none','FaceAlpha',opac);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[73, 135, 59]/255,'LineStyle','none','FaceAlpha',opac);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[169, 208, 88]/255,'LineStyle','none','FaceAlpha',opac);
f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[203, 143, 70]/255,'LineStyle','none','FaceAlpha',opac);

axis equal
axis off
%view([a1 b1])
view([[11.1385,-68.9693]])
camlight('left')
material dull 
hold on
%scatter3(newX1(l0,1) + rand(length(l0),1)*1,newY1(l0,1)+ rand(length(l0),1)*1,newZ1(l0,1)+ rand(length(l0),1)*1,W1(l0)*60,[12,156,245]/255,'filled')
scatter3(newX1(l1,1) + rand(length(l1),1)*5,newY1(l1,1)+ rand(length(l1),1)*5,newZ1(l1,1)+ rand(length(l1),1)*5,W1(l1)*60,[45, 122, 234]/255,'filled')
scatter3(newX1(l2,1) + rand(length(l2),1)*5,newY1(l2,1)+ rand(length(l2),1)*5,newZ1(l2,1)+ rand(length(l2),1)*5,W1(l2)*60,[33, 105, 192]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
scatter3(newX1(l4,1) + rand(length(l4),1)*5,newY1(l4,1)+ rand(length(l4),1)*5,newZ1(l4,1)+ rand(length(l4),1)*5,W1(l4)*60,[255, 140, 100]/255,'filled')
scatter3(newX1(l5,1) + rand(length(l5),1)*5,newY1(l5,1)+ rand(length(l5),1)*5,newZ1(l5,1)+ rand(length(l5),1)*5,W1(l5)*60,[46, 137, 46]/255,'filled')

%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[117,76,36]/255,'filled')
scatter3(newX1(l6,1) + rand(length(l6),1)*1,newY1(l6,1)+ rand(length(l6),1)*1,newZ1(l6,1)+ rand(length(l6),1)*1,W1(l6)*60,[214, 140, 52]/255,'filled')
scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,[226, 193, 32]/255,'filled')
scatter3(newX1(l8,1) + rand(length(l8),1)*1,newY1(l8,1)+ rand(length(l8),1)*1,newZ1(l8,1)+ rand(length(l8),1)*1,W1(l8)*60,[226, 193, 32]/255,'filled')


colormap(parula);
print(h,['~/Desktop/Am_CS5_AllTissues_withXiang2.png'],'-dpng','-r600')









colormap(parula);
print(h,['~/Desktop/Am_CS5_AllTissues_withXiang3.png'],'-dpng','-r600')









close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
%f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l1,1) + rand(length(l1),1)*1,newY1(l1,1)+ rand(length(l1),1)*1,newZ1(l1,1)+ rand(length(l1),1)*1,W1(l1)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
colormap(parula);
print(h,['~/Desktop/Am_CS6_Embissues.png'],'-dpng','-r600')


%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
colormap(parula);
print(h,['~/Desktop/AmEmDisc_CS6_AllTissues.png'],'-dpng','-r600')

close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
%f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
%C1 = readtable(File2);
%CAno = C1(l2,1)
%textscatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,table2array(CAno),'TextDensityPercentage',100)

colormap(parula);
print(h,['~/Desktop/AmEmDisc_CS6_Embissues.png'],'-dpng','-r600')





%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
colormap(parula);
print(h,['~/Desktop/EmDisc_CS6_AllTissues.png'],'-dpng','-r600')

close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
%f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
%C1 = readtable(File2);
%CAno = C1(l3,1)
%textscatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,table2array(CAno),'TextDensityPercentage',100)
colormap(parula);
print(h,['~/Desktop/EmDisc_CS6_Embissues.png'],'-dpng','-r600')

close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',1);
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
%f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
colormap(parula);
print(h,['~/Desktop/TissueRef3.png'],'-dpng','-r600')








%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l4,1) + rand(length(l4),1)*1,newY1(l4,1)+ rand(length(l4),1)*1,newZ1(l4,1)+ rand(length(l4),1)*1,W1(l4)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
colormap(parula);
print(h,['~/Desktop/ExMes_CS6_AllTissues.png'],'-dpng','-r600')

close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
f8 = patch('Faces',OBJ2.objects(36).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);

axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l4,1) + rand(length(l4),1)*1,newY1(l4,1)+ rand(length(l4),1)*1,newZ1(l4,1)+ rand(length(l4),1)*1,W1(l4)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
colormap(parula);
print(h,['~/Desktop/ExMes_CS6_Embissues.png'],'-dpng','-r600')









%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l5,1) + rand(length(l5),1)*1,newY1(l5,1)+ rand(length(l5),1)*1,newZ1(l5,1)+ rand(length(l5),1)*1,W1(l5)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
colormap(parula);
print(h,['~/Desktop/ExMesSYS_CS6_AllTissues.png'],'-dpng','-r600')





%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l6,1) + rand(length(l6),1)*1,newY1(l6,1)+ rand(length(l6),1)*1,newZ1(l6,1)+ rand(length(l6),1)*1,W1(l6)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
colormap(parula);
print(h,['~/Desktop/VE_CS6_AllTissues.png'],'-dpng','-r600')

close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
%f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l6,1) + rand(length(l6),1)*1,newY1(l6,1)+ rand(length(l6),1)*1,newZ1(l6,1)+ rand(length(l6),1)*1,W1(l6)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
%C1 = readtable(File2);
%CAno = C1(l6,1)
%textscatter3(newX1(l6,1) + rand(length(l6),1)*1,newY1(l6,1)+ rand(length(l6),1)*1,newZ1(l6,1)+ rand(length(l6),1)*1,table2array(CAno))
colormap(parula);
print(h,['~/Desktop/VE_CS6_Embissues.png'],'-dpng','-r600')


close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',1);
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
%f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
colormap(parula);
print(h,['~/Desktop/TissueRef2.png'],'-dpng','-r600')



%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
colormap(parula);
print(h,['~/Desktop/PGC_CS6_AllTissues.png'],'-dpng','-r600')

close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,[230, 230, 0]/255,'filled')
scatter3(newX1(l9,1) + rand(length(l9),1)*1,newY1(l9,1)+ rand(length(l9),1)*1,newZ1(l9,1)+ rand(length(l9),1)*1,W1(l9)*60,[230, 230, 0]/255,'filled')
scatter3(newX1(l10,1) + rand(length(l10),1)*1,newY1(l10,1)+ rand(length(l10),1)*1,newZ1(l10,1)+ rand(length(l10),1)*1,W1(l10)*60,[230, 230, 0]/255,'filled')


scatter3(newX1(l1,1) + rand(length(l1),1)*1,newY1(l1,1)+ rand(length(l1),1)*1,newZ1(l1,1)+ rand(length(l1),1)*1,W1(l1)*60,[155, 17, 171]/255,'filled')
scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[103, 99, 212]/255,'filled')
scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[22, 17, 171]/255,'filled')

%scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l9,1) + rand(length(l9),1)*1,newY1(l9,1)+ rand(length(l9),1)*1,newZ1(l9,1)+ rand(length(l9),1)*1,W1(l9)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l10,1) + rand(length(l10),1)*1,newY1(l10,1)+ rand(length(l10),1)*1,newZ1(l10,1)+ rand(length(l10),1)*1,W1(l10)*60,[96,56,19]/255,'filled')
%scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,[12,156,245]/255,'filled')
%C1 = readtable(File2);
%CAno = C1(l7,1)
%textscatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,table2array(CAno))
%C1 = readtable(File2);
%CAno = C1(l10,1)
%textscatter3(newX1(l10,1) + rand(length(l10),1)*1,newY1(l10,1)+ rand(length(l10),1)*1,newZ1(l10,1)+ rand(length(l10),1)*1,table2array(CAno),'TextDensityPercentage',100)
%colormap(parula);


print(h,['~/Desktop/PGC_CS6_Embissues_ANO.png'],'-dpng','-r600')




close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,[12,156,245]/255,'filled')
scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[12,156,245]/255,'filled')
scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')

%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
colormap(parula);
print(h,['~/Desktop/PGC_CS6_Embissues.png'],'-dpng','-r600')


close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
colormap(parula);
print(h,['~/Desktop/PGC_TissueRef.png'],'-dpng','-r600')




%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l8,1) + rand(length(l8),1)*1,newY1(l8,1)+ rand(length(l8),1)*1,newZ1(l8,1)+ rand(length(l8),1)*1,W1(l8)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
colormap(parula);
print(h,['~/Desktop/SYS_CS6_AllTissues.png'],'-dpng','-r600')

close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
%f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l8,1) + rand(length(l8),1)*1,newY1(l8,1)+ rand(length(l8),1)*1,newZ1(l8,1)+ rand(length(l8),1)*1,W1(l8)*60,[12,156,245]/255,'filled')
%scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[240,76,4]/255,'filled')
%scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[96,56,19]/255,'filled')
%C1 = readtable(File2);
%CAno = C1(l8,1)
%textscatter3(newX1(l8,1) + rand(length(l8),1)*1,newY1(l8,1)+ rand(length(l8),1)*1,newZ1(l8,1)+ rand(length(l8),1)*1,table2array(CAno))

colormap(parula);
print(h,['~/Desktop/SYS_CS6_Embissues.png'],'-dpng','-r600')








close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[60,221,221]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[45,126,224]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[225,84,232]/255,'LineStyle','none','FaceAlpha',.1);
%f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[255,140,100]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[46,137,46]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on
scatter3(newX1(l8,1) + rand(length(l8),1)*1,newY1(l8,1)+ rand(length(l8),1)*1,newZ1(l8,1)+ rand(length(l8),1)*1,W1(l8)*60,[6,137,46]/255,'filled')
scatter3(newX1(l1,1) + rand(length(l1),1)*1,newY1(l1,1)+ rand(length(l1),1)*1,newZ1(l1,1)+ rand(length(l1),1)*1,W1(l1)*60,[60,221,221]/255,'filled')
scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[45,126,224]/255,'filled')
scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[45,126,224]/255,'filled')
scatter3(newX1(l6,1) + rand(length(l6),1)*1,newY1(l6,1)+ rand(length(l6),1)*1,newZ1(l6,1)+ rand(length(l6),1)*1,W1(l6)*60,[255,140,100]/255,'filled')
scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,[225,84,232]/255,'filled')

colormap(parula);
print(h,['~/Desktop/SYSVEEMDISC_CS6_Embissuesv2.png'],'-dpng','-r600')





close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[60,221,221]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[45,126,224]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[225,84,232]/255,'LineStyle','none','FaceAlpha',.1);
%f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[255,140,100]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[46,137,46]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
%scatter3(newX1(find(isnan(newX1)==0),1),newY1(find(isnan(newX1)==0),1),newZ1(find(isnan(newX1)==0),1),W1(find(isnan(newX1)==0))*60,'k','filled')
hold on

%SYS,Am,Am/EmDisc,EmDisc,VE

scatter3(newX1(l8,1) + rand(length(l8),1)*1,newY1(l8,1)+ rand(length(l8),1)*1,newZ1(l8,1)+ rand(length(l8),1)*1,W1(l8)*60,[6,137,46]/255,'filled')
scatter3(newX1(l1,1) + rand(length(l1),1)*1,newY1(l1,1)+ rand(length(l1),1)*1,newZ1(l1,1)+ rand(length(l1),1)*1,W1(l1)*60,[60,221,221]/255,'filled')
scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[45,126,224]/255,'filled')
scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[45,126,224]/255,'filled')
scatter3(newX1(l6,1) + rand(length(l6),1)*1,newY1(l6,1)+ rand(length(l6),1)*1,newZ1(l6,1)+ rand(length(l6),1)*1,W1(l6)*60,[255,140,100]/255,'filled')

scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,[225,84,232]/255,'filled')


colormap(parula); 
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gcf, 'InvertHardcopy', 'off')
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print(h,['~/Desktop/SYSVEEMDISC_CS6_Embissues_black.png'],'-dpng','-r600')





f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[60,221,221]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[45,126,224]/255,'LineStyle','none','FaceAlpha',1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[225,84,232]/255,'LineStyle','none','FaceAlpha',1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[159,209,63]/255,'LineStyle','none','FaceAlpha',1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[255,140,100]/255,'LineStyle','none','FaceAlpha',1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[46,137,46]/255,'LineStyle','none','FaceAlpha',1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);




h=figure(1)

scatter(1,1,100,[179, 123, 34]/255,'filled')
hold on
scatter(1,2,100,[155, 17, 171]/255,'filled')
scatter(1,3,100,[103, 99, 212]/255,'filled')
scatter(1,4,100,[22, 17, 171]/255,'filled')
scatter(1,5,100,[181, 36, 51]/255,'filled')
colormap(parula);
ylim([0 2])
ylim([0 6])
print(h,['~/Desktop/Colours.svg'],'-svg','-r600')




close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[60,221,221]/255,'LineStyle','none','FaceAlpha',1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[45,126,224]/255,'LineStyle','none','FaceAlpha',1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[225,84,232]/255,'LineStyle','none','FaceAlpha',1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[159,209,63]/255,'LineStyle','none','FaceAlpha',1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[255,140,100]/255,'LineStyle','none','FaceAlpha',1);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[46,137,46]/255,'LineStyle','none','FaceAlpha',1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
colormap(parula);
print(h,['~/Desktop/TissueRef_noTb.png'],'-dpng','-r600')



close all
%EmDisc
h=figure(1)
%subplot(2,4,1);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[60,221,221]/255,'LineStyle','none','FaceAlpha',0);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[45,126,224]/255,'LineStyle','none','FaceAlpha',1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[225,84,232]/255,'LineStyle','none','FaceAlpha',0);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[159,209,63]/255,'LineStyle','none','FaceAlpha',0);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[255,140,100]/255,'LineStyle','none','FaceAlpha',0);
f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[46,137,46]/255,'LineStyle','none','FaceAlpha',0);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([36.2185,-67.4732])
view([29.3426,-90])
camlight('left')
material dull 
hold on
colormap(parula);
print(h,['~/Desktop/TissueRef_noTb_emb.png'],'-dpng','-r600')



[OBJ1,section] = LoadCS5('3D');
[DCS5,LocationsCS5,XYZCS5,CellTypeCS5,ShotsCS5] = LoadShots('CS5');

File1 = '/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/S1_set1_all.csv'
File2 = '/Users/christopherpenfold/Desktop/Data/Endometrial/InVitro/Matteo/C1_set1.csv'
%[newX1,newY1,newZ1,X1,Y1,Z1,W1,CellType1,CellTypeReg] = ProjectData(File1,File2,ShotsCS6,10);

[newX1,newY1,newZ1,X1,Y1,Z1,W1,CellType1] = ProjectDataCS5(File1,File2,ShotsCS5,10);
% l1 = length(newX1(find(isnan(newX1)==0),1));
%
%l1 = find( strcmp(CellType1,'EmDisc_d14')==1 | strcmp(CellType1,'Am_d14')==1  |  strcmp(CellType1,'Am/EmDisc_d14')==1  | strcmp(CellType1,'Am/EmDisc')==1 );
%l1a = find( strcmp(CellType1,'EmDisc_d14')==1  );
%l1b = find( strcmp(CellType1,'Am_d14')==1  );
%l1c = find(  strcmp(CellType1,'Am/EmDisc_d14')==1  );
%l1d = find(  strcmp(CellType1,'Am/EmDisc')==1 );
%l3 = find( strcmp(CellType1,'ExMes_d14')==1 | strcmp(CellType1,'putPGC')==1 | strcmp(CellType1,'putExMes')==1 |  strcmp(CellType1,'Hyp/Am')==1);
%l2 = find( strcmp(CellType1,'Hyp_d14')==1 );

writecell(CellType1,'~/Desktop/CT_CS5.csv') 
writecell(CellTypeReg,'~/Desktop/CTR_CS5.csv') 

l1 = find( isnan(newX1)==0 & (strcmp(CellType1,'Am')==1 ) );
l2 = find( isnan(newX1)==0 & (strcmp(CellType1,'Am/EmDisc')==1 ) );
l3 = find( isnan(newX1)==0 & (strcmp(CellType1,'EmDisc')==1 ) );
l4 = find( isnan(newX1)==0 & (strcmp(CellType1,'ExMes')==1 ) );
l5 = find( isnan(newX1)==0 & (strcmp(CellType1,'ExMes/SYS')==1 ) );
l6 = find( isnan(newX1)==0 & (strcmp(CellType1,'Hypoblast')==1 ) );

l7 = find( isnan(newX1)==0 & (strcmp(CellType1,'PGC')==1 ) );
l8 = find( isnan(newX1)==0 & (strcmp(CellType1,'SYS')==1 ) );


h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([285.5574,-0.4695])
camlight('right')
hold on
scatter3(newX1(l1,1) + rand(length(l1),1)*1,newY1(l1,1)+ rand(length(l1),1)*1,newZ1(l1,1)+ rand(length(l1),1)*1,W1(l1)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/EmAm1_CS5_Embissues.png'],'-dpng','-r600')

close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',.1);
%%f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,56,19]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([3.6634,-86.9955])
camlight('right')
hold on
scatter3(newX1(l1,1) + rand(length(l1),1)*1,newY1(l1,1)+ rand(length(l1),1)*1,newZ1(l1,1)+ rand(length(l1),1)*1,W1(l1)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/EmAm1_CS5.png'],'-dpng','-r600')




close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([285.5574,-0.4695])
camlight('right')
hold on
scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/AmEmDisc_CS5_Embissues.png'],'-dpng','-r600')

close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',.1);
%%f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,56,19]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([3.6634,-86.9955])
camlight('right')
hold on
scatter3(newX1(l2,1) + rand(length(l2),1)*1,newY1(l2,1)+ rand(length(l2),1)*1,newZ1(l2,1)+ rand(length(l2),1)*1,W1(l2)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/AmEmDisc_CS5.png'],'-dpng','-r600')







close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([285.5574,-0.4695])
camlight('right')
hold on
scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/EmDisc_CS5_Embissues.png'],'-dpng','-r600')

close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',.1);
%%f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,56,19]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([3.6634,-86.9955])
camlight('right')
hold on
scatter3(newX1(l3,1) + rand(length(l3),1)*1,newY1(l3,1)+ rand(length(l3),1)*1,newZ1(l3,1)+ rand(length(l3),1)*1,W1(l3)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/EmDisc_CS5.png'],'-dpng','-r600')






close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([285.5574,-0.4695])
camlight('right')
hold on
scatter3(newX1(l4,1) + rand(length(l4),1)*1,newY1(l4,1)+ rand(length(l4),1)*1,newZ1(l4,1)+ rand(length(l4),1)*1,W1(l4)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/ExMes_CS5_Embissues.png'],'-dpng','-r600')

close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',.1);
%%f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,56,19]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([3.6634,-86.9955])
camlight('right')
hold on
scatter3(newX1(l4,1) + rand(length(l4),1)*1,newY1(l4,1)+ rand(length(l4),1)*1,newZ1(l4,1)+ rand(length(l4),1)*1,W1(l4)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/ExMes_CS5.png'],'-dpng','-r600')





close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([285.5574,-0.4695])
camlight('right')
hold on
scatter3(newX1(l5,1) + rand(length(l5),1)*1,newY1(l5,1)+ rand(length(l5),1)*1,newZ1(l5,1)+ rand(length(l5),1)*1,W1(l5)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/ExMesSYS_CS5_Embissues.png'],'-dpng','-r600')

close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',.1);
%%f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,56,19]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([3.6634,-86.9955])
camlight('right')
hold on
scatter3(newX1(l5,1) + rand(length(l5),1)*1,newY1(l5,1)+ rand(length(l5),1)*1,newZ1(l5,1)+ rand(length(l5),1)*1,W1(l5)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/ExMesSYS_CS5.png'],'-dpng','-r600')






close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,56,19]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([285.5574,-0.4695])
camlight('right')
hold on
scatter3(newX1(l6,1) + rand(length(l6),1)*1,newY1(l6,1)+ rand(length(l6),1)*1,newZ1(l6,1)+ rand(length(l6),1)*1,W1(l6)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/Hypoblast_CS5_Embissues.png'],'-dpng','-r600')

close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',.1);
%%f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,56,19]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([3.6634,-86.9955])
camlight('right')
hold on
scatter3(newX1(l6,1) + rand(length(l6),1)*1,newY1(l6,1)+ rand(length(l6),1)*1,newZ1(l6,1)+ rand(length(l6),1)*1,W1(l6)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/Hypoblast_CS5.png'],'-dpng','-r600')






close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([285.5574,-0.4695])
camlight('right')
hold on
scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,'r','filled')
scatter3(newX1(l9,1) + rand(length(l9),1)*1,newY1(l9,1)+ rand(length(l9),1)*1,newZ1(l9,1)+ rand(length(l9),1)*1,W1(l9)*60,'g','filled')
scatter3(newX1(l10,1) + rand(length(l10),1)*1,newY1(l10,1)+ rand(length(l10),1)*1,newZ1(l10,1)+ rand(length(l10),1)*1,W1(l10)*60,'b','filled')

print(h,['~/Desktop/PGC_CS5_Embissues.png'],'-dpng','-r600')

close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',.1);
%%f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,56,19]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([3.6634,-86.9955])
camlight('right')
hold on
scatter3(newX1(l7,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,W1(l7)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/PGC_CS5.png'],'-dpng','-r600')






close all
h=figure(2)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);

axis equal
axis off
material dull 
colormap(parula);
view([285.5574,-0.4695])
camlight('right')
hold on
scatter3(newX1(l8,1) + rand(length(l8),1)*1,newY1(l8,1)+ rand(length(l8),1)*1,newZ1(l8,1)+ rand(length(l8),1)*1,W1(l8)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/SYS_CS5_Embissues.png'],'-dpng','-r600')

close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',.1);
%%f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,56,19]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([3.6634,-86.9955])
camlight('right')
hold on
scatter3(newX1(l8,1) + rand(length(l8),1)*5,newY1(l8,1)+ rand(length(l8),1)*5,newZ1(l8,1)+ rand(length(l8),1)*5,W1(l8)*600,[12,156,245]/255,'filled')
print(h,['~/Desktop/SYS_CS5.png'],'-dpng','-r600')

figure(1)
C1 = readtable(File2);
CAno = C1(l8,1)
textscatter3(newX1(l8,1) + rand(length(l7),1)*1,newY1(l7,1)+ rand(length(l7),1)*1,newZ1(l7,1)+ rand(length(l7),1)*1,table2array(CAno))




close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',.1);

axis equal
axis off
material dull 
colormap(parula);
view([285.5574,-0.4695])
camlight('right')
hold on
scatter3(newX1(l11,1) + rand(length(l11),1)*1,newY1(l11,1)+ rand(length(l11),1)*1,newZ1(l11,1)+ rand(length(l11),1)*1,W1(l11)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/SYS_CS5_Embissues.png'],'-dpng','-r600')

close all
h=figure(1)
f1 = patch('Faces',OBJ1.objects(4).data.vertices,'Vertices',  OBJ1.vertices,'FaceColor',[135,123,214]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ1.objects(20).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[12, 156, 245]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ1.objects(8).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[240,76,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ1.objects(12).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[230,134,0]/255,'LineStyle','none','FaceAlpha',.1);
%%f5 = patch('Faces',OBJ1.objects(24).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,31,230]/255,'LineStyle','none','FaceAlpha',.1);
f6 = patch('Faces',OBJ1.objects(16).data.vertices,'Vertices',OBJ1.vertices,'FaceColor',[96,56,19]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
material dull 
colormap(parula);
view([3.6634,-86.9955])
camlight('right')
hold on
scatter3(newX1(l11,1) + rand(length(l11),1)*1,newY1(l11,1)+ rand(length(l11),1)*1,newZ1(l11,1)+ rand(length(l11),1)*1,W1(l11)*60,[12,156,245]/255,'filled')
print(h,['~/Desktop/SYS_CS5.png'],'-dpng','-r600')






[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,ShotsCS6] = LoadShots('CS6');

File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set1.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set1.csv'
[newX1,newY1,newZ1,X1,Y1,Z1,W1,CellType1,RegCellType1] = ProjectData(File1,File2,ShotsCS6,5);

[newX1c,newY1c,newZ1c] = CorrectProjectData(newX1,newY1,newZ1,W1,CellType1,RegCellType1,OBJ2);

File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set2.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set2.csv'
[newX2,newY2,newZ2,X2,Y2,Z2,W2,CellType2] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set3.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set3.csv'
[newX3,newY3,newZ3,X3,Y3,Z3,W3,CellType3] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set4.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set4.csv'
[newX4,newY4,newZ4,X4,Y4,Z4,W4,CellType4] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set5.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set5.csv'
[newX5,newY5,newZ5,X5,Y5,Z5,W5,CellType5] = ProjectData(File1,File2,ShotsCS6,5);

File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set1.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set1.csv'
[newXs1,newYs1,newZs1,Xs1,Ys1,Zs1,Ws1,CellTypes1] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set2.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set2.csv'
[newXs2,newYs2,newZs2,Xs2,Ys2,Zs2,Ws2,CellTypes2] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set3.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set3.csv'
[newXs3,newYs3,newZs3,Xs3,Ys3,Zs3,Ws3,CellTypes3] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set4.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set4.csv'
[newXs4,newYs4,newZs4,Xs4,Ys4,Zs4,Ws4,CellTypes4] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set5.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set5.csv'
[newXs5,newYs5,newZs5,Xs5,Ys5,Zs5,Ws5,CellTypes5] = ProjectData(File1,File2,ShotsCS6,5);

S1 = sqrt( (newX1-X1).^2 + (newY1-Y1).^2 + (newZ1-Z1).^2 );
S2 = sqrt( (newX2-X2).^2 + (newY2-Y2).^2 + (newZ2-Z2).^2 );
S3 = sqrt( (newX3-X3).^2 + (newY3-Y3).^2 + (newZ3-Z3).^2 );
S4 = sqrt( (newX4-X4).^2 + (newY4-Y4).^2 + (newZ4-Z4).^2 );
S5 = sqrt( (newX5-X5).^2 + (newY5-Y5).^2 + (newZ5-Z5).^2 );

S1_ = sqrt( (newXs1-Xs1).^2 + (newYs1-Ys1).^2 + (newZs1-Zs1).^2 );
S2_ = sqrt( (newXs2-Xs2).^2 + (newYs2-Ys2).^2 + (newZs2-Zs2).^2 );
S3_ = sqrt( (newXs3-Xs3).^2 + (newYs3-Ys3).^2 + (newZs3-Zs3).^2 );
S4_ = sqrt( (newXs4-Xs4).^2 + (newYs4-Ys4).^2 + (newZs4-Zs4).^2 );
S5_ = sqrt( (newXs5-Xs5).^2 + (newYs5-Ys5).^2 + (newZs5-Zs5).^2 );

RME = [S1,S1_,S2,S2_,S3,S3_,S4,S4_,S5,S5_];

csvwrite('MSE_N=5.csv',RME)
T = cell2table([CellType1,CellType2,CellType3,CellType4,CellType5,CellTypes1,CellTypes2,CellTypes3,CellTypes4,CellTypes5])
writetable(T,'CellType_N=5.csv')


ind1_1 = find(strcmp(CellType1,'EmDisc')==1);
ind1_2 = find(strcmp(CellType2,'EmDisc')==1);
ind1_3 = find(strcmp(CellType3,'EmDisc')==1);
ind1_4 = find(strcmp(CellType4,'EmDisc')==1);
ind1_5 = find(strcmp(CellType5,'EmDisc')==1);

inds1_1 = find(strcmp(CellTypes1,'EmDisc')==1);
inds1_2 = find(strcmp(CellTypes2,'EmDisc')==1);
inds1_3 = find(strcmp(CellTypes3,'EmDisc')==1);
inds1_4 = find(strcmp(CellTypes4,'EmDisc')==1);
inds1_5 = find(strcmp(CellTypes5,'EmDisc')==1);

ind2_1 = find(strcmp(CellType1,'VE')==1);
ind2_2 = find(strcmp(CellType2,'VE')==1);
ind2_3 = find(strcmp(CellType3,'VE')==1);
ind2_4 = find(strcmp(CellType4,'VE')==1);
ind2_5 = find(strcmp(CellType5,'VE')==1);

inds2_1 = find(strcmp(CellTypes1,'VE')==1);
inds2_2 = find(strcmp(CellTypes2,'VE')==1);
inds2_3 = find(strcmp(CellTypes3,'VE')==1);
inds2_4 = find(strcmp(CellTypes4,'VE')==1);
inds2_5 = find(strcmp(CellTypes5,'VE')==1);


ind3_1 = find(strcmp(CellType1,'Am')==1);
ind3_2 = find(strcmp(CellType2,'Am')==1);
ind3_3 = find(strcmp(CellType3,'Am')==1);
ind3_4 = find(strcmp(CellType4,'Am')==1);
ind3_5 = find(strcmp(CellType5,'Am')==1);

inds3_1 = find(strcmp(CellTypes1,'Am')==1);
inds3_2 = find(strcmp(CellTypes2,'Am')==1);
inds3_3 = find(strcmp(CellTypes3,'Am')==1);
inds3_4 = find(strcmp(CellTypes4,'Am')==1);
inds3_5 = find(strcmp(CellTypes5,'Am')==1);


ind4_1 = find(strcmp(CellType1,'Tb')==1);
ind4_2 = find(strcmp(CellType2,'Tb')==1);
ind4_3 = find(strcmp(CellType3,'Tb')==1);
ind4_4 = find(strcmp(CellType4,'Tb')==1);
ind4_5 = find(strcmp(CellType5,'Tb')==1);

inds4_1 = find(strcmp(CellTypes1,'Tb')==1);
inds4_2 = find(strcmp(CellTypes2,'Tb')==1);
inds4_3 = find(strcmp(CellTypes3,'Tb')==1);
inds4_4 = find(strcmp(CellTypes4,'Tb')==1);
inds4_5 = find(strcmp(CellTypes5,'Tb')==1);

%EmDisc
h=figure(1)
subplot(2,4,1);
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])
camlight('left')
material dull 
hold on
scatter3(newX1(ind1_1,1),newY1(ind1_1,1),newZ1(ind1_1,1),W1(ind1_1)*60,'k','filled')
scatter3(newX2(ind1_2,1),newY1(ind1_2,1),newZ1(ind1_2,1),W2(ind1_2)*60,'k','filled')
scatter3(newX3(ind1_3,1),newY1(ind1_3,1),newZ1(ind1_3,1),W3(ind1_3)*60,'k','filled')
scatter3(newX4(ind1_4,1),newY1(ind1_4,1),newZ1(ind1_4,1),W4(ind1_4)*60,'k','filled')
scatter3(newX5(ind1_5,1),newY1(ind1_5,1),newZ1(ind1_5,1),W5(ind1_5)*60,'k','filled')
%scatter3(newXs1(inds1_1,1),newYs1(inds1_1,1),newZs1(inds1_1,1),Ws1(inds1_1)*60,'r','filled')
%scatter3(newXs2(inds1_2,1),newYs1(inds1_2,1),newZs1(inds1_2,1),Ws2(inds1_2)*60,'r','filled')
%scatter3(newXs3(inds1_3,1),newYs1(inds1_3,1),newZs1(inds1_3,1),Ws3(inds1_3)*60,'r','filled')
%scatter3(newXs4(inds1_4,1),newYs1(inds1_4,1),newZs1(inds1_4,1),Ws4(inds1_4)*60,'r','filled')
%scatter3(newXs5(inds1_5,1),newYs1(inds1_5,1),newZs1(inds1_5,1),Ws5(inds1_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);

subplot(2,4,5);
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newXs1(inds1_1,1),newYs1(inds1_1,1),newZs1(inds1_1,1),Ws1(inds1_1)*60,'r','filled')
scatter3(newXs2(inds1_2,1),newYs1(inds1_2,1),newZs1(inds1_2,1),Ws2(inds1_2)*60,'r','filled')
scatter3(newXs3(inds1_3,1),newYs1(inds1_3,1),newZs1(inds1_3,1),Ws3(inds1_3)*60,'r','filled')
scatter3(newXs4(inds1_4,1),newYs1(inds1_4,1),newZs1(inds1_4,1),Ws4(inds1_4)*60,'r','filled')
scatter3(newXs5(inds1_5,1),newYs1(inds1_5,1),newZs1(inds1_5,1),Ws5(inds1_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);


subplot(2,4,2);
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newX1(ind2_1,1),newY1(ind2_1,1),newZ1(ind2_1,1),W1(ind2_1)*60,'k','filled')
scatter3(newX2(ind2_2,1),newY1(ind2_2,1),newZ1(ind2_2,1),W2(ind2_2)*60,'k','filled')
scatter3(newX3(ind2_3,1),newY1(ind2_3,1),newZ1(ind2_3,1),W3(ind2_3)*60,'k','filled')
scatter3(newX4(ind2_4,1),newY1(ind2_4,1),newZ1(ind2_4,1),W4(ind2_4)*60,'k','filled')
scatter3(newX5(ind2_5,1),newY1(ind2_5,1),newZ1(ind2_5,1),W5(ind2_5)*60,'k','filled')
%scatter3(newXs1(inds1_1,1),newYs1(inds1_1,1),newZs1(inds1_1,1),Ws1(inds1_1)*60,'r','filled')
%scatter3(newXs2(inds1_2,1),newYs1(inds1_2,1),newZs1(inds1_2,1),Ws2(inds1_2)*60,'r','filled')
%scatter3(newXs3(inds1_3,1),newYs1(inds1_3,1),newZs1(inds1_3,1),Ws3(inds1_3)*60,'r','filled')
%scatter3(newXs4(inds1_4,1),newYs1(inds1_4,1),newZs1(inds1_4,1),Ws4(inds1_4)*60,'r','filled')
%scatter3(newXs5(inds1_5,1),newYs1(inds1_5,1),newZs1(inds1_5,1),Ws5(inds1_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);

subplot(2,4,6);
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newXs1(inds2_1,1),newYs1(inds2_1,1),newZs1(inds2_1,1),Ws1(inds2_1)*60,'r','filled')
scatter3(newXs2(inds2_2,1),newYs1(inds2_2,1),newZs1(inds2_2,1),Ws2(inds2_2)*60,'r','filled')
scatter3(newXs3(inds2_3,1),newYs1(inds2_3,1),newZs1(inds2_3,1),Ws3(inds2_3)*60,'r','filled')
scatter3(newXs4(inds2_4,1),newYs1(inds2_4,1),newZs1(inds2_4,1),Ws4(inds2_4)*60,'r','filled')
scatter3(newXs5(inds2_5,1),newYs1(inds2_5,1),newZs1(inds2_5,1),Ws5(inds2_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);








subplot(2,4,3);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newX1(ind3_1,1),newY1(ind3_1,1),newZ1(ind3_1,1),W1(ind3_1)*60,'k','filled')
scatter3(newX2(ind3_2,1),newY1(ind3_2,1),newZ1(ind3_2,1),W2(ind3_2)*60,'k','filled')
scatter3(newX3(ind3_3,1),newY1(ind3_3,1),newZ1(ind3_3,1),W3(ind3_3)*60,'k','filled')
scatter3(newX4(ind3_4,1),newY1(ind3_4,1),newZ1(ind3_4,1),W4(ind3_4)*60,'k','filled')
scatter3(newX5(ind3_5,1),newY1(ind3_5,1),newZ1(ind3_5,1),W5(ind3_5)*60,'k','filled')
%scatter3(newXs1(inds1_1,1),newYs1(inds1_1,1),newZs1(inds1_1,1),Ws1(inds1_1)*60,'r','filled')
%scatter3(newXs2(inds1_2,1),newYs1(inds1_2,1),newZs1(inds1_2,1),Ws2(inds1_2)*60,'r','filled')
%scatter3(newXs3(inds1_3,1),newYs1(inds1_3,1),newZs1(inds1_3,1),Ws3(inds1_3)*60,'r','filled')
%scatter3(newXs4(inds1_4,1),newYs1(inds1_4,1),newZs1(inds1_4,1),Ws4(inds1_4)*60,'r','filled')
%scatter3(newXs5(inds1_5,1),newYs1(inds1_5,1),newZs1(inds1_5,1),Ws5(inds1_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);

subplot(2,4,7);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newXs1(inds3_1,1),newYs1(inds3_1,1),newZs1(inds3_1,1),Ws1(inds3_1)*60,'r','filled')
scatter3(newXs2(inds3_2,1),newYs1(inds3_2,1),newZs1(inds3_2,1),Ws2(inds3_2)*60,'r','filled')
scatter3(newXs3(inds3_3,1),newYs1(inds3_3,1),newZs1(inds3_3,1),Ws3(inds3_3)*60,'r','filled')
scatter3(newXs4(inds3_4,1),newYs1(inds3_4,1),newZs1(inds3_4,1),Ws4(inds3_4)*60,'r','filled')
scatter3(newXs5(inds3_5,1),newYs1(inds3_5,1),newZs1(inds3_5,1),Ws5(inds3_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);

subplot(2,4,4);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
%view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newX1(ind4_1,1),newY1(ind4_1,1),newZ1(ind4_1,1),W1(ind4_1)*60,'k','filled')
scatter3(newX2(ind4_2,1),newY1(ind4_2,1),newZ1(ind4_2,1),W2(ind4_2)*60,'k','filled')
scatter3(newX3(ind4_3,1),newY1(ind4_3,1),newZ1(ind4_3,1),W3(ind4_3)*60,'k','filled')
scatter3(newX4(ind4_4,1),newY1(ind4_4,1),newZ1(ind4_4,1),W4(ind4_4)*60,'k','filled')
scatter3(newX5(ind4_5,1),newY1(ind4_5,1),newZ1(ind4_5,1),W5(ind4_5)*60,'k','filled')
%scatter3(newXs1(inds1_1,1),newYs1(inds1_1,1),newZs1(inds1_1,1),Ws1(inds1_1)*60,'r','filled')
%scatter3(newXs2(inds1_2,1),newYs1(inds1_2,1),newZs1(inds1_2,1),Ws2(inds1_2)*60,'r','filled')
%scatter3(newXs3(inds1_3,1),newYs1(inds1_3,1),newZs1(inds1_3,1),Ws3(inds1_3)*60,'r','filled')
%scatter3(newXs4(inds1_4,1),newYs1(inds1_4,1),newZs1(inds1_4,1),Ws4(inds1_4)*60,'r','filled')
%scatter3(newXs5(inds1_5,1),newYs1(inds1_5,1),newZs1(inds1_5,1),Ws5(inds1_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);

subplot(2,4,8);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newXs1(inds4_1,1),newYs1(inds4_1,1),newZs1(inds4_1,1),Ws1(inds4_1)*60,'r','filled')
scatter3(newXs2(inds4_2,1),newYs1(inds4_2,1),newZs1(inds4_2,1),Ws2(inds4_2)*60,'r','filled')
scatter3(newXs3(inds4_3,1),newYs1(inds4_3,1),newZs1(inds4_3,1),Ws3(inds4_3)*60,'r','filled')
scatter3(newXs4(inds4_4,1),newYs1(inds4_4,1),newZs1(inds4_4,1),Ws4(inds4_4)*60,'r','filled')
scatter3(newXs5(inds4_5,1),newYs1(inds4_5,1),newZs1(inds4_5,1),Ws5(inds4_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);


%EmDisc
%h=figure(1)
subplot(2,4,2);
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
%view([c1,d1])
camlight('left')
material dull 
hold on
for i = 1:length(ind1)
%plot3([newX(ind1(i),1);X(ind1(i),1)],[newY(ind1(i),1);X(ind1(i),1)],[newZ(ind1(i),1);Z(ind1(i),1)],'k.-')
scatter3(X(ind1(i),1),Y(ind1(i),1),Z(ind1(i),1),W(ind1(i))*60,'k','filled')
end
for i = 1:length(ind8)
%plot3([newX(ind1(i),1);X(ind1(i),1)],[newY(ind1(i),1);X(ind1(i),1)],[newZ(ind1(i),1);Z(ind1(i),1)],'k.-')
scatter3(X(ind8(i),1),Y(ind8(i),1),Z(ind8(i),1),W(ind8(i))*60,'r','filled')
end
for i = 1:length(ind7)
%plot3([newX(ind1(i),1);X(ind1(i),1)],[newY(ind1(i),1);X(ind1(i),1)],[newZ(ind1(i),1);Z(ind1(i),1)],'k.-')
scatter3(X(ind7(i),1),Y(ind7(i),1),Z(ind7(i),1),W(ind7(i))*60,'b','filled')
end
colormap(parula);


%Density mapping
allemdisc=[newX1(ind1_1,1),newY1(ind1_1,1),newZ1(ind1_1,1);
newX2(ind1_2,1),newY1(ind1_2,1),newZ1(ind1_2,1);
newX3(ind1_3,1),newY1(ind1_3,1),newZ1(ind1_3,1);
newX4(ind1_4,1),newY1(ind1_4,1),newZ1(ind1_4,1);
newX5(ind1_5,1),newY1(ind1_5,1),newZ1(ind1_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to EmDisc')
ax = gca
nClim = ax.CLim;
ax.CLim = nClim;

subplot(2,4,2);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,3);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to Am')
ax = gca
ax.CLim = nClim;

subplot(2,4,4);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to Tb')
ax = gca
ax.CLim = nClim;


allemdisc=[newXs1(ind1_1,1),newYs1(ind1_1,1),newZs1(ind1_1,1);
newXs2(ind1_2,1),newYs1(ind1_2,1),newZs1(ind1_2,1);
newXs3(ind1_3,1),newYs1(ind1_3,1),newZs1(ind1_3,1);
newXs4(ind1_4,1),newYs1(ind1_4,1),newZs1(ind1_4,1);
newXs5(ind1_5,1),newYs1(ind1_5,1),newZs1(ind1_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,5);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to EmDisc')
ax = gca
ax.CLim = nClim;

subplot(2,4,6);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,7);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to Am')
ax = gca
ax.CLim = nClim;

subplot(2,4,8);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to Tb')
ax = gca
ax.CLim = nClim;

print(h,['~/Desktop/EmDisc.png'],'-dpng','-r600')
close all


%Density mapping
allemdisc=[newX1(ind2_1,1),newY1(ind2_1,1),newZ1(ind2_1,1);
newX2(ind2_2,1),newY1(ind2_2,1),newZ1(ind2_2,1);
newX3(ind2_3,1),newY1(ind2_3,1),newZ1(ind2_3,1);
newX4(ind2_4,1),newY1(ind2_4,1),newZ1(ind2_4,1);
newX5(ind2_5,1),newY1(ind2_5,1),newZ1(ind2_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to EmDisc')
ax = gca
nClim = ax.CLim;
ax.CLim = nClim;

subplot(2,4,2);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,3);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to Am')
ax = gca
ax.CLim = nClim;


subplot(2,4,4);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to Tb')
ax = gca
ax.CLim = nClim;



allemdisc=[newXs1(ind2_1,1),newYs1(ind2_1,1),newZs1(ind2_1,1);
newXs2(ind2_2,1),newYs1(ind2_2,1),newZs1(ind2_2,1);
newXs3(ind2_3,1),newYs1(ind2_3,1),newZs1(ind2_3,1);
newXs4(ind2_4,1),newYs1(ind2_4,1),newZs1(ind2_4,1);
newXs5(ind2_5,1),newYs1(ind2_5,1),newZs1(ind2_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,5);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to EmDisc')
ax = gca
ax.CLim = nClim;

subplot(2,4,6);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,7);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to Am')
ax = gca
ax.CLim = nClim;

subplot(2,4,8);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to Tb')
ax = gca
ax.CLim = nClim;


%Now Amnion
%Density mapping
allemdisc=[newX1(ind3_1,1),newY1(ind3_1,1),newZ1(ind3_1,1);
newX2(ind3_2,1),newY1(ind3_2,1),newZ1(ind3_2,1);
newX3(ind3_3,1),newY1(ind3_3,1),newZ1(ind3_3,1);
newX4(ind3_4,1),newY1(ind3_4,1),newZ1(ind3_4,1);
newX5(ind3_5,1),newY1(ind3_5,1),newZ1(ind3_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to EmDisc')
ax = gca
nClim = ax.CLim;
ax.CLim = nClim;

subplot(2,4,2);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,3);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to Am')
ax = gca
ax.CLim = nClim;


subplot(2,4,4);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to Tb')
ax = gca
ax.CLim = nClim;



allemdisc=[newXs1(ind3_1,1),newYs1(ind3_1,1),newZs1(ind3_1,1);
newXs2(ind3_2,1),newYs1(ind3_2,1),newZs1(ind3_2,1);
newXs3(ind3_3,1),newYs1(ind3_3,1),newZs1(ind3_3,1);
newXs4(ind3_4,1),newYs1(ind3_4,1),newZs1(ind3_4,1);
newXs5(ind3_5,1),newYs1(ind3_5,1),newZs1(ind3_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,5);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to EmDisc')
ax = gca
ax.CLim = nClim;

subplot(2,4,6);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,7);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to Am')
ax = gca
ax.CLim = nClim;

subplot(2,4,8);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to Tb')
ax = gca
ax.CLim = nClim;




%Now Tb
%Density mapping
allemdisc=[newX1(ind4_1,1),newY1(ind4_1,1),newZ1(ind4_1,1);
newX2(ind4_2,1),newY1(ind4_2,1),newZ1(ind4_2,1);
newX3(ind4_3,1),newY1(ind4_3,1),newZ1(ind4_3,1);
newX4(ind4_4,1),newY1(ind4_4,1),newZ1(ind4_4,1);
newX5(ind4_5,1),newY1(ind4_5,1),newZ1(ind4_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to EmDisc')
ax = gca
nClim = ax.CLim;
ax.CLim = nClim;

subplot(2,4,2);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,3);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to Am')
ax = gca
ax.CLim = nClim;


subplot(2,4,4);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to Tb')
ax = gca
ax.CLim = nClim;



allemdisc=[newXs1(ind4_1,1),newYs1(ind4_1,1),newZs1(ind4_1,1);
newXs2(ind4_2,1),newYs1(ind4_2,1),newZs1(ind4_2,1);
newXs3(ind4_3,1),newYs1(ind4_3,1),newZs1(ind4_3,1);
newXs4(ind4_4,1),newYs1(ind4_4,1),newZs1(ind4_4,1);
newXs5(ind4_5,1),newYs1(ind4_5,1),newZs1(ind4_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,5);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to EmDisc')
ax = gca
ax.CLim = nClim;

subplot(2,4,6);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,7);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to Am')
ax = gca
ax.CLim = nClim;

subplot(2,4,8);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to Tb')
ax = gca
ax.CLim = nClim;
