
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

parser <- OptionParser(
    usage       = "Rscript %prog output",
    description = "\nan merge lagged ipmr simulations",
    epilogue    = ""
  )


cli <- parse_args(parser, positional_arguments = 1)

output <- cli$args[1]

df <- as.list(file.path(output, list.files(output)))
df <- lapply(df, readRDS) %>% bind_rows

str(df)

strings <- stringr::str_split(list.files(output)[1], "[[:punct:]]")
strings 

out.file = paste(output,strings[[1]][3],"_", strings[[1]][4], "_", "Merged.RDS", sep = "")
out.file 

saveRDS(df, file = out.file)
