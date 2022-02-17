function set_vbr_path(vbr_path)
    fileID = fopen('./vbr_path.txt','w');
    fprintf(fileID,'%s', vbr_path);
    fclose(fileID);
end
