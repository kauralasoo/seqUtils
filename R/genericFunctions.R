idVectorToList <- function(id_vector){
  #Convert a vector of IDs into a list of ids where the elements have the same name
  id_list = as.list(id_vector)
  names(id_list) = id_vector
  return(id_list)
}