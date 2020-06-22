
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


stringr::str_split(list.files(output)[1], "[[:punct:]]")

saveRDS(df, file = file.path(output, "Merged.RDS"))