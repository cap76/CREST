addpath(genpath('/Users/christopherpenfold/Desktop/Code/Spatial\ modelling/SpatialModelling/'))
addpath(genpath('/Users/christopherpenfold/Desktop/Code/Spatial modelling/SpatialModelling'))
cd('/Users/christopherpenfold/Desktop/Code/Spatial modelling/SpatialModelling')

%Load CS6
[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,ShotsCS6] = LoadShots('CS6');

File1 = '/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping5/S1_set1_all_withPreImpl_byCl.csv'
File2 = '/Users/christopherpenfold/Desktop/Data/MatteoMarmosetMapping5/C1_set1_withPreImpl_byCl.csv'
[newX1,newY1,newZ1,X1,Y1,Z1,W1,CellType1,CellTypeReg] = ProjectData(File1,File2,ShotsCS6,10);

l0 = find( isnan(newX1)==0 );
S1 = readtable(File1);
C1 = readtable(File2);

l1 = find( isnan(newX1)==0 & (strcmp(CellType1,'Amnion')==1 ) );
l2 = find( isnan(newX1)==0 & (strcmp(CellType1,'EmbDisc')==1 ) );
l3 = find( isnan(newX1)==0 & (strcmp(CellType1,'ExEmbMes')==1 ) );
l4 = find( isnan(newX1)==0 & (strcmp(CellType1,'Hypoblast')==1 ) );
l5 = find( isnan(newX1)==0 & (strcmp(CellType1,'YolkSac')==1 ) );
l6 = find( isnan(newX1)==0 & (strcmp(CellType1,'STB')==1 ) );
l7 = find( isnan(newX1)==0 & (strcmp(CellType1,'CTB')==1 ) );
l8 = find( isnan(newX1)==0 & (strcmp(CellType1,'eCTB')==1 ) );

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
print(h,['~/Desktop/PGC_CS6_Embissues_ANO.png'],'-dpng','-r600')




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
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'Color','k')
set(gcf, 'color', 'k')
set(gcf, 'InvertHardcopy', 'off')
fig=gcf;ax=fig.CurrentAxes;fig.Color='k';fig.OuterPosition=fig.InnerPosition;
axis off
print(h,['~/Desktop/SYSVEEMDISC_CS6_Embissues_black.png'],'-dpng','-r600')


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
