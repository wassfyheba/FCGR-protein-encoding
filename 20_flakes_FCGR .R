
library (seqinr)

prot=read.fasta('protein.fasta',seqtype = "AA",as.string = FALSE)
prot2=unlist(prot)
library(kaos)

prot2.cgr=cgr(prot2,res=100)
prot3=cgr.plot(prot2.cgr,mode='points')
prot4=cgr.plot(prot2.cgr,mode='matrix')
prot5=cgr.res(prot4, 200)
prot6=vectorize(prot2.cgr)






