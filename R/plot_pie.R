#' Plot pie
#'
#' @param vector A numberic vector for plot
#' @param expression A Strings for calculate percent
#' @param label A character vector for display catalog
#'
#' @return NULL
#' @export
#'
plot_pie <- function(vector,expression,label){
  per <- c(mean(eval(parse(text = expression))),1-mean(eval(parse(text = expression))))
  pie(per,labels = paste0(label," (",round(per,digits = 3)*100,"%)"),
           border="white", col=c("#66C2A5", "#FC8D62"))
  gridGraphics::grid.echo()
  a <- grid.grab()
  return(a)
}
