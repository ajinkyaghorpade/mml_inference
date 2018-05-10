prepare_data_mlogit <- function(mydata)
{
  
  # Create pooled panel data from cross sectional data input
  panel.mlogit.data <- format_long_data(mydata);
  
  return(mydata);
}


format_wide_data <- function(data) {
  panel.y_dummy <- t(sapply(lapply(data$y, function(x) x[[1]]), unlist));
  panel.y <- max.col(panel.y_dummy);
  panel.x <- t(sapply(lapply(data$X, function(x) x[[1]]), unlist));
  panel.data.frame <- as.data.frame(cbind(panel.x,panel.y));
  names(panel.data.frame)[1:(dim(data$X[[1]][[1]])[2] * dim(data$X[[1]][[1]])[1])] <- paste(sort(rep(seq(1,dim(data$X[[1]][[1]])[2]), dim(data$X[[1]][[1]])[1])),rep(seq(1,dim(data$X[[1]][[1]])[2]), dim(data$X[[1]][[1]])[1]), sep = ".");
  panel.mlogit.data <- mlogit.data(panel.data.frame, choice = "panel.y", shape = "wide", varying = 1:9);
  return (panel.mlogit.data);
}

format_long_data <- function(data) {
  panel.x <- do.call(rbind, lapply(data$X,function(x) do.call(rbind, x)));
  panel.x <- lapply(data$X, function(t) lapply(t, function(t) as.data.frame(t)));
  panel.x <- do.call(rbind, lapply(panel.x,function(x) rbindlist(x, idcol = 'id')));
  panel.x <- rbindlist(lapply(panel.x,function(x) rbindlist(x, idcol = 'id')), idcol = 'ind.id');
  panel.x <- rbindlist(lapply(data$X,function(x) rbindlist(x,idcol = 'id')));
  panel.y <- do.call(c, lapply(data$y,function(x) do.call(c, x)));
  panel.alt_idx <- do.call(c, lapply(data$y, function(x) do.call(c, lapply(x, function(t) do.call(c, list(seq(1,length(t))))))));
  panel.data.frame <- as.data.frame(cbind(panel.x,panel.y, panel.alt_idx));
  panel.mlogit.data <- mlogit.data(panel.data.frame, choice = "panel.y", shape = "long", alt.var = "panel.alt_idx");
  return(panel.mlogit.data);
}

# TODO: Try faster rearrangement of the lists
# library(data.table)
# head(rbindlist(data))