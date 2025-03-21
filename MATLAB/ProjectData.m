function [newX,newY,newZ,X,Y,Z,W,CellType,CellTypeReg] = ProjectData(File1,File2,ShotsCS6,N)
%This is the options.


S1 = readtable(File1);
C1 = readtable(File2);


targs = C1(1:end,2) %.Idents_mammal_combined__inds1_ %textdata(2:end,2);
regs = C1(1:end,14:23);


targs = strrep(targs.Variables,'P1_','')
targs = strrep(targs,'P2_','')

regs = strrep([regs(:,1).Variables,regs(:,2).Variables,regs(:,3).Variables,regs(:,4).Variables,regs(:,5).Variables,regs(:,6).Variables,regs(:,7).Variables,regs(:,8).Variables,regs(:,9).Variables,regs(:,10).Variables],'P1_','')
%regs = strrep(regs,'P1_','')
regs = strrep(regs,'P2_','')

for i = 1:length( ShotsCS6.textdata(2:end,1) )
     ShotsCS6.textdata{i+1,1} = ['E15C2_' ShotsCS6.textdata{i+1,1}];
end

targInd = zeros(length(targs),1);

X = NaN*ones(length(targs),1);
Y = NaN*ones(length(targs),1);
Z = NaN*ones(length(targs),1);

CellType = cell(length(targs),1);

for i = 1:length(targs)
    try
        
      %  if i==170
       %     keyboard
        targInd(i,1) = find(strcmp(targs{i},  ShotsCS6.textdata(2:end,1) )==1);        
        X(i,1) = ShotsCS6.data(targInd(i,1),1);
        Y(i,1) = ShotsCS6.data(targInd(i,1),2);
        Z(i,1) = ShotsCS6.data(targInd(i,1),3);              
        CellType{i,1} = ShotsCS6.textdata{ targInd(i,1)+1 , 2};
        %keyboard
        end
    %end
end

%%%%%Added in for when not comparing to our own SS2 dataset
%CellType = C1.textdata(2:end,3);
CellType = C1(1:end,3);
CellType = CellType.Variables




regInd = zeros(size(regs,1),size(regs,2));

Xt = NaN*ones( size(regs,1),size(regs,2) );
Yt = NaN*ones( size(regs,1),size(regs,2) );
Zt = NaN*ones( size(regs,1),size(regs,2) );

CellTypeReg = cell(size(regs,1),size(regs,2));


for i = 1:size(regs,1)
    %if i==170
    %    keyboard
    for j = 1:size(regs,2)
        try
            
        regInd(i,j) = find(strcmp(regs{i,j},  ShotsCS6.textdata(2:end,1) )==1);
        Xt(i,j) = ShotsCS6.data(regInd(i,j),1);
        Yt(i,j) = ShotsCS6.data(regInd(i,j),2);
        Zt(i,j) = ShotsCS6.data(regInd(i,j),3);   
         CellTypeReg{i,j} = ShotsCS6.textdata{ regInd(i,j)+1 , 2};
       
        end
    end
    %end
end

W = NaN*ones(length(targs),1);
newX = NaN*ones(length(targs),1);
newY = NaN*ones(length(targs),1);
newZ = NaN*ones(length(targs),1);

Dat = [S1.Var2,S1.Var3,S1.Var4,S1.Var5,S1.Var6,S1.Var7,S1.Var8,S1.Var9,S1.Var10,S1.Var11];

for i = 1:size(X,1)
    
   % if i==170
   %     keyboard
   u = Dat(i,1:N);
   iX = Xt(i,1:N);
    
    inds = find(isnan(iX)==0) % & u~=0);
    if isempty(inds)==0
        v = u(inds);
        v = v./sum(v);
        newX(i,1) = sum(v.*Xt(i,inds));
        newY(i,1) = sum(v.*Yt(i,inds));
        newZ(i,1) = sum(v.*Zt(i,inds));                                
        W(i,1) = sum(u(inds))/sum(u);        
    end
    %end
end