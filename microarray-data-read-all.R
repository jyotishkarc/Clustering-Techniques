
library(dplyr)
library(stringr)


path <- "D:/My Documents/Datasets/CompCancer Database/TwoClass/"

files <- list.files(path)

for(h in 1:length(files)){
   
   dataset.name <- files[h] %>% str_remove("_database.xlsx") %>% str_replace_all("-","_")
   
   dataset <- paste0(path, files[h]) %>% 
      read_excel() %>% 
      as.matrix() %>%
      apply(c(1,2), function(val) as.numeric(val))
   
   assign(paste0(dataset.name,".gt"), dataset[ ,1] %>% unlist())
   assign(paste0(dataset.name,".data"), dataset[ ,-1])
   
   print(h)
}

