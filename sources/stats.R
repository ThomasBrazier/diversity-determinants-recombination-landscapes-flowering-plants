# Absolute Cross-validation
# Lee & Cox 2010
ACV = function(data = data.frame()) {
  criterion = (1/nrow(data))*sum(abs(data$obs - data$pred), na.rm = TRUE)
  return(criterion)
}

# Least squares Cross validation
# Lee & Cox 2010
LSCV = function(data = data.frame()) {
  criterion = (1/nrow(data))*sum((data$obs - data$pred)^2, na.rm = TRUE)
  return(criterion)
}