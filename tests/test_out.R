df<-data.frame(Column1 = c("1","2","3"), Column2 = c("a","b","c"))
colnames(df)<-c("num", "chr")

write.csv(df,"/home/users6/sshukor/scripts/test_1.csv", row.names = FALSE)
write.csv(df,"\\home\\users6\\sshukor\\scripts\\test_2.csv", row.names = FALSE)