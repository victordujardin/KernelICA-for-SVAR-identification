## Find negative impacts

neg.ica       <- vector("list",length(index.ica))
for (i in 1:length(index.ica)){
  neg.ica[[i]]       <- which(as.vector(A.ID.ica[[index.ica[i]]]) <0)
}
neg.kica       <- vector("list",length(index.kica))
for (i in 1:length(index.ica)){
  neg.kica[[i]]       <- which(as.vector(A.ID.kica[[index.kica[i]]]) <0)
}

neg.dc        <- vector("list",length(index.dc))
for (i in 1:length(index.dc)){
  neg.dc[[i]]       <- which(as.vector(A.ID.dc[[index.dc[i]]]) <0)
}
neg.pml       <- vector("list",length(index.pml))
for (i in 1:length(index.pml)){
  neg.pml[[i]]       <- which(as.vector(A.ID.pml[[index.pml[i]]]) <0)
}

