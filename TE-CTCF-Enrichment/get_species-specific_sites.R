df1 = read.table("tmp1", header=FALSE, sep="\t", stringsAsFactors=FALSE)
df2 = read.table("tmp2", header=FALSE, sep="\t", stringsAsFactors=FALSE)
rownames(df1) = df1$V4
rownames(df2) = df2$V4
df3 = df1[which(!(rownames(df1) %in% rownames(df2))),]
df3$V5 = 0
df3$V6 = 1
write.table(df3, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, file="tmp3")
