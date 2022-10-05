

#==========================================================
# Saving logs
#==========================================================
# All operations on genetic and Marey maps must be saved in a log file, in a single line with with date and parameters, for reproducibility purposes
# path: "", path where to save the log file, which will take exactly the same name as the genetic/Marey map with .log extension, in a log directory
# if "", save in a log directory in the exact path of the Marey map
# msg: "", the message to be written, if "" then the call of the function is saved
# append: TRUE, passed to write.file, set FALSE if you want to reinitialise the file
save.log = function(path = "/data/Marey_maps", msg = "", append = TRUE) {
  # Check if 'log' directory exists, otherwise create the dir
  if (!dir.exists(paste(wd, path, "/log", sep = ""))) {
    dir.create(paste(wd, path, "/log", sep = ""))
  }
  
  # If msg = "", then insert the call to the initial function
  if (msg == "") {
      msg = deparse(sys.calls()[[sys.nframe()-1]])
  }
  # Add a time (day + hours:minutes) to the log text
  time = date()
  # The log phrase is a new line appended to a '.log' file in a specific log directory
  log.msg = paste("[", set, " : ", time, "] ", msg, sep = "")
  #print(log.msg)
  write.table(log.msg, file = paste(wd, path, "/log/", gsub(".txt", "", set), ".log", sep = ""), append = append,
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

resetlog = function(path = "/data/Marey_maps", set = "") {
  file.remove(paste(path, set, sep = ""))
}

#==========================================================
# Plot Marey maps
#==========================================================
# Plot a Marey map dataframe with ggplot and custom theme
plot.Marey = function() {
  
}



#==========================================================
# latin2arabic
#==========================================================
# Convert latin numerals to arabic numbers in a character vector

# x is a vector of character type

latin2arabic = function(x) {
  
  replacelatin = function(char) {
    if (grepl("chromosome[ ]*[IVX]+", char)) {
      # Find chromosome name in numerals
      replacement = regmatches(as.character(char), gregexpr("chromosome [IVX]*", as.character(char)))
      replacement = gsub("chromosome ", "", gsub(",", "", replacement))
      # Convert roman to numeric
      replacement =  utils:::.roman2numeric(as.character(replacement))
      char = gsub("chromosome [IVX]*", paste("chromosome ", replacement, sep = ""), char)
    }
    return(char)
  }
  # Apply to the whole vector and return to user
  results = unlist(lapply(x, replacelatin))
  return(results)
}
