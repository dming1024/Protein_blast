
#获取氨基酸序列
#1. 将gene symbol 转换为 Uniprotkb ID
#2. 根据uniprokb ID，和API进行检索
# curl -X GET "https://rest.uniprot.org/uniprotkb/O15305" -H "accept: application/json"
#3. 从检索结果中获取序列信息

gs=c("CNOT6L","CALM2","CCNI","SPINT2")

library('org.Hs.eg.db')
annotation <- select(org.Hs.eg.db, 
                          keys=gs, 
                          columns=c('UNIPROT', 'SYMBOL', 'ENTREZID'), keytype="SYMBOL")

install.packages("UniprotR",dependencies = T)
library(UniprotR)
Accessions <-GetAccessionList("https://s3.amazonaws.com/csvpastebin/uploads/9571fa356c67a0c7c95e8431799a051a/Accessions.csv") 
TaxaObj <- GetNamesTaxa(Accessions) 
TaxaObj
PlotChromosomeInfo(TaxaObj)

ConvertID(Accessions,ID_from = "ACC+ID" , ID_to = "Gene Name")

GETSeqFastaUniprot(Accessions[1],FilePath = "./", FileName = "abc.fa")
