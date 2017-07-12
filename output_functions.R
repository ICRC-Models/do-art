## Function that takes single vector and converts it to data.frame using the original array as a reference.

reconstructArray <- function(vec) {
  array(data = vec,
        dim = c(length(age), length(sex), length(hiv)),
        dimnames = list(age_names, sex_names, hiv_names))
  
}

array2DF <- function(a, df_names) {
  
  require(plyr)
  df <- adply(a, c(1:3), function(x) {
    
    y <- data.frame(cbind(x, rownames(x)))
    return(y)
  })
  
  names(df) <- df_names
  
  return(df)
  
}

test <- reconstruct(pop)
array2DF(test, df_names = c("age", "sex", "hiv", "n"))

z <- apply(out[1:4,], 1, function(x) {
  
  ## print(x)
##   print(class(x))
  
  arr <- reconstructArray(x[2:length(x)])
  ## print(arr)
  df <- array2DF(arr, df_names = c("age", "sex", "hiv", "n"))
  df$time <- x[1]
  
  return(df)
  
  
  
})

tt <- adply(test, c(1:3), function(x) {
  
  y <- data.frame(cbind(x, rownames(x)))
  return(y)
  
})

names(tt) <- c("age", "sex", "hiv", "n")

a <- expand.grid(age = age_names, sex = sex_names, hiv = hiv_names)
a$n <- c(pop)
