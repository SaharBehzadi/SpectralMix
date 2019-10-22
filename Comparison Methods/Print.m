function Print(filepath,algorithm,result)
filename=strcat(filepath,'\',algorithm,'label.txt');
fileID = fopen(filename,'w');
fprintf(fileID, '%i \n', result);
fclose(fileID);