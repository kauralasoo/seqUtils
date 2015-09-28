saveEqtlbmaData <- function(data_list, output_dir, file_suffix = "exression", 
                            project_root = "/nfs/users/nfs_k/ka8/group-scratch/kaur/projects/macrophage-gxe-study/"){
  #Save data for eqtlbma to disk
  
  #Construct file list
  file_list = data_frame(sample_name = names(data_list), 
                         file_path = file.path(project_root, output_dir, 
                                               paste(names(data_list),file_suffix, "txt.gz", sep = ".")))
  file_list_path = file.path(output_dir, 
                             paste("list",file_suffix, "txt", sep = "."))
  write.table(file_list, file_list_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  #Save each matrix as a separate gzipped txt file
  print(file_list)
  for (sn in file_list$sample_name){
    file_path = file.path(output_dir, paste(sn,file_suffix, "txt.gz", sep = "."))
    print(file_path)
    file = gzfile(file_path, "w")
    write.table(data_list[[sn]], file, sep = "\t", quote = FALSE)
    close(file)
  }
}