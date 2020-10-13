
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

parser <- OptionParser(
    usage       = "Rscript %prog output",
    description = "\nan merge lagged ipmr simulations",
    epilogue    = ""
  )


cli <- parse_args(parser, positional_arguments = 1)

output <- cli$args[1]

file <- as.list(file.path(output, list.files(output)))

output
str(file)

df <- lapply(file, function(x) readRDS(x) %>% 
               select(-contains("hist"))) %>% 
      bind_rows

#df1 <- lapply(file, function(x) readRDS(x) %>%
#              select(contains("hist"))) %>%
#      bind_rows

str(df)
#str(df1)

strings <- stringr::str_split(list.files(output)[1], "[[:punct:]]")
strings 

out.file = paste(output,strings[[1]][3],"_", strings[[1]][4], "_", "Merged.RDS", sep = "")
#out.file1 = paste(output,strings[[1]][3],"_", strings[[1]][4], "_", "Hist.RDS", sep = "")

out.file 
#out.file1

saveRDS(df, file = out.file)
#saveRDS(df1, file = out.file1)
