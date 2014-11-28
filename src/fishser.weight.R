fisher.score <- function(Col, y)
{
  uni <-  unique(y)
  Col1 <-  Col[y==uni[1]]
  Col2 <-  Col[y==uni[2]]
  num <- length(Col1) * (mean(Col1) - mean(Col) ) +  length(Col2) * (mean(Col2) - mean(Col) ) 
  den <-  length(Col1) * var(Col1) +   length(Col2) * var(Col2)
 num/den 
}
