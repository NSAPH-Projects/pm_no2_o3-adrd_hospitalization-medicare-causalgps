#
# Objective: A template processing file for the project.
# In this template, we want to have access to external data and plot some of
# them.
#
# Original author: Naeem Khoshnevis
# Contributing authors: N/A
# Last update: February 28, 2023

# Notes:
# Each sub projects is a R or Rmd file that has three sections:
#     1 - Initialization, which requires sub project name and whether you want
#         to symlink external_public and external_private or not.
#     2 - Body, this is the place that you put your code. Use path_obj to have
#         access to sub project related path.
#     3 - Finalizing, which documents the used packages + hash of the used code.

# initiate the process ---------------------------------------------------------
file_name_env <- sys.frame(1)
sp_name  <- sub("\\.R$", "", basename(file_name_env$fileName))
path_obj <- initialize_sub_project(sp_name = sp_name)

# setup cache on disk
cdb <- cachem::cache_disk(path_obj$sp_cache)

# ------------------------------------------------------------------------------

# Processing -------------------------------------------------------------------
# Write tmp function to load data
load_data_tmp <- function(data_path){
  # read in full data; 34,763,397 rows
  data <- fst::read_fst(data_path,
                        as.data.table = TRUE)
  data.table::setnames(data,
                       old = c("pct_blk", "pct_owner_occ"),
                       new = c("prop_blk", "prop_owner_occ"))
  data[, `:=`(zip = as.factor(zip),
              year = as.factor(year),
              cohort = as.factor(cohort),
              age_grp = as.factor(age_grp),
              sex = as.factor(sex),
              race = as.factor(race),
              dual = as.factor(dual))]

  # get small sub sample
  data <- data[1:100000,]

  return(data)
}

# memoise function
m_load_data_tmp <- memoise::memoise(load_data_tmp, cache = cdb)

# load data
data_path <- file.path(path_obj$dir_data_private_ext_1,
                       "main_data_1_private",
                       "analysis",
                       "ADRD_complete_tv.fst")

data <- m_load_data_tmp(data_path = data_path)


pdf(file.path(path_obj$sp_output,"density_plot.pdf"), width = 8, height = 8)
plot(density(data$tmmx))
dev.off()

# Finialize --------------------------------------------------------------------
finalize_subproject(path_obj = path_obj)
