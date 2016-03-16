## conversion char to num
c2num <- function(x){
    as.numeric(gsub(",", ".", x, fixed = TRUE))
}
