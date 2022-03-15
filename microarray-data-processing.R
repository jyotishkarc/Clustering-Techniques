dataset <- khan.2001

temp <- cbind(t(dataset)[-1,-c(1:2)], labels.rename(t(dataset)[-1,2]))
temp2 <- apply(temp, c(1,2), function(val) as.numeric(val))

khan.2001 <- temp2

write_xlsx(as.data.frame(temp2), 
           path = "D:\\My Documents\\microarray-data\\khan-2001.xlsx")
