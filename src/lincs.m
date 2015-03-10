path('/cbcbhomes/vinash85/software/l1ktools/matlab/', path);
path('/cbcbhomes/vinash85/software/l1ktools/matlab/lib', path);
lincs = parse_gctx('~/project/project2-scratch/data/lincs/l1000/level2/gex_epsilon_n1429794x978.gctx')
// tab=tdfread('~/project/project2-scratch/data/lincs/l1000/metadata/inst.info');
// save('~/project/project2-scratch/data/lincs/l1000/level2/gex_deltap_n53094x978.mat', lincs)
lincsMat = lincs.mat;
rid=lincs.rid;
cid=lincs.cid; 
// save('~/project/project2-scratch/data/lincs/l1000/level2/gex_epsilon_n1429794x978.dim.mat', 'cid', 'rid');
save('~/project/project2-scratch/data/lincs/l1000/level2/gex_epsilon_n1429794x978.cid.mat', 'cid');
save('~/project/project2-scratch/data/lincs/l1000/level2/gex_epsilon_n1429794x978.rid.mat', 'rid')
// csvwrite('~/project/project2-scratch/data/lincs/l1000/level2/gex_epsilon_n1429794x978.matrix.csv', lincsMat);
save('~/project/project2-scratch/data/lincs/l1000/level2/gex_epsilon_n1429794x978.matrix.mat', 'lincsMat');




###writing to file 
[nrows,ncols] = size(lincs.cid);
fileID = fopen('lincs.cid', 'w');
formatSpec = '%s\n';
for row = 1:nrows
    fprintf(fileID,formatSpec,lincs.cid{row});
end

[nrows,ncols] = size(lincs.rid);
fileID = fopen('lincs.rid', 'w');
formatSpec = '%s\n';
for row = 1:nrows
    fprintf(fileID,formatSpec,lincs.rid{row});
end

csvwrite('lincs.sklMat.csv', lincMat);