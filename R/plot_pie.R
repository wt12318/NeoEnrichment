#' Plot pie
#'
#' Plot pie plot for a numeric vector
#' @param vector A numberic vector for plot
#' @param expression A String, can be parsed as expression for calculate percent, usually a conditional expression
#' @param label A character vector , for display catalog
#'
#' @return NULL
#' @export
#' @examples
#' a <- c(1,2,3,4)
#' plot_pie(a,expression = "a < 2",c("a < 2","a >= 2"))
#'
plot_pie <- function(vector,expression,label){
  per <- c(mean(eval(parse(text = expression))),1-mean(eval(parse(text = expression))))
  a <- pie(per,labels = paste0(label," (",round(per,digits = 3)*100,"%)"),
           border="white", col=c("#66C2A5", "#FC8D62"))
  return(a)
}
