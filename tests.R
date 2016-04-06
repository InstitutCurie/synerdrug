d <- read.table("../longData.txt", header = TRUE)

object <- new("DrugSyn", data = d, doses = list(A=unique(d$A), B=unique(d$B)))
