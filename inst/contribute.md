# How to Contribute
_Revision 1.0.0_        
_Last update: March 1, 2023_


We follow the idea of compendium proposed by _Robert Gentleman & Duncan Temple Lang_ in **Statistical Analyses and Reproducible Research** paper. 

## Objective

When designing this workflow and structure the following objectives were kept in mind:

- It should be possible for all research contributors to easily access any commit within the project and run the code without requiring extensive modifications.
- Each user's results should be stored in their own account or project folder. In the lab account, working on a shared project folder facilitates collaboration on ongoing research; however, it increases the chance of writing results in others dedicated project folders. 
- Since we experiment with various hypotheses in our research, it is important to make setting up a run and collecting results as straightforward as possible without writing too much code by the researchers.
- Each run should be documented, and the documentation should include details of the environment for each sub-project, as well as the code used to carry out the research.
- Access to data beyond the project folder should be seamless. 
- Users should be able to access as many external data sets as they wish from different resources.

## Strucutre

The project has a similar R package structure. However, only some requirements may be necessary to apply. We start with minimum requirements and will update this page for more upcoming conventions.

The following shows an example of the structure of the project, files, and folders.

```sh
.
├── ADRDMedicare220121.Rproj
├── study_data
│   ├── private
│   │   ├── external -> path/to/external/data
│   │   └── internal
│   └── public
│       └── internal
├── DESCRIPTION
├── inst
│   └── contribute.md
├── man
│   ├── document_code_version.Rd
│   ├── finalize_subproject.Rd
│   ├── initialize_sub_project.Rd
│   ├── setup_path.Rd
│   └── setup_subproj_folders.Rd
├── NAMESPACE
├── NEWS.md
├── R
│   ├── document_code_version.R
│   ├── finalize_subproject.R
│   ├── initialize_sub_project.R
│   ├── setup_path.R
│   └── setup_subproj_folders.R
├── README.md
├── results
│   └── inp_2302_plot_curated_data
│       ├── cache_files
│       │   ├── ac491ee9ad9017279cbf182bf8b50693.rds
│       │   ├── c619389318be0ed22f3b8ce43a25d230.rds
│       │   └── c76b8384592081a913a69113c4bab59a.rds
│       ├── log_files
│       └── output_files
│           ├── density_plot.pdf
│           └── sessioninfo_20230228190525.txt
└── rscripts
    └── inp_230228_template.R

```
In the following, we discuss how to work with and contribute to the package. 

## Setting up working environment

- If you are working with pure R (e.g., with batch submissions), you need to install the package and run the functions. [TODO: This is not tested]
- In case of working with RStudio, follow these steps:
  - We use `devtools` to load, test, and generate documentation. 
  - Make sure you are in the project folder, or click on `ADRDMedicare220121.Rproj` to activate the project.
  - Run `devtools::load_all()`
  - If you have added a new function or modified the documentation, run `devtools::document()`
- Each user should have a `R/external_path` file in the project root directory. This is personalized information and should not be committed to git. In the following, we will discuss how to use it. 
  - If any unit tests have been added to the project, run `devtools::test()` and make sure that it does not fail. Any function that is being added under `R` folder should have unit tests to make sure that some other member or yourself does not modify and break it. 
  - Run `devtools::check()` to make sure that the package is in good shape and does not violate the general package structure.

## Running a scientific workflow

After setting up the environment and making sure that
  - `devtools::test()`, and 
  - `devtools::load_all()`
are running without a problem. We can start conducting research. All `.R` and `.Rmd` files go into the `rscript` folder. After implementing all functions and getting a feeling that the function can be used in multiple places or by other users, you can add more functions under `R` folder. Please consider adding unit tests for those functions or at least add them to TODO list.

Before doing any analyses, remove the `external_[number]` folders in the `study_data/private` and the `study_data/public` folders. You may get warnings otherwis. [TODO: write a code to do this internally. Make sure that recursive soft link does not cause problem at the target.]

We call each scientific task a sub-project. A sub-project can be any number of steps that takes some input values and data and provide some outputs. Each sub-project, which is either `.R` or `.Rmd` files, has the following structure (all in the same file). [TODO: the `.Rmd` files have not been tested yet.]

### Objective of sub-project

One paragraph about the objective of the subproject needs to be added. This should include enough information to help yourself and others to understand why you are doing this processing. If you need more details, open an issue, add more details, and put the issue number in this subproject section. In the case of an `.R` file, this part will be a comment.

```r
#  [First section general info]
# 
#  Your description about the subproject goes here.
#
# Original Author Name:  Your name
# Contributors Name: Contributors name [Contributor names of just this file, not the entire project.]
# Last update: MM-DD-YYYY  
# 
# Also mention if the project needs external data or not. You do not need to provide full path, however, you can mention the data, who curated it.
#
```

### Sub-project name

The second section is just one line of code, which is the sub-project name. The prefix and date after the convention help others to find the sub-projects and sort them out quickly. It is recommended to choose the same file name and project name. 

The name of the sub-project starts with three character prefix, then date (YYMMDD), and then any descriptive name. A folder with the same name will be created under the `results` folder, and all related documents will be there.

The convention for prefixes:
  - inp: in progress (any exploratory code)
  - man: manuscript (any code that directly or indirectly generates data or figures for manuscript)
  - or pick any other prefix and update the convention.

```r
# [Second section]
sp_name <- "inp_230228_[any name]"

```

### Sub-project initialization

Each sub-project requires initialization. The initialization sets the appropriate path internally, so you can rest assured that all required information will be stored in dedicated folders. Use the following command to initialize the sub-project.

```r
# [Third section]
path_obj <- initialize_sub_project(sp_name = sp_name)
```

`sp_name` is the sub-project name that you defined in the previous step. Please note that all access to data in the code follows the following structure. 

```sh
├── study_data
│   ├── private
│   │   ├── external_1 -> path/to/external/data
│   │   └── internal
│   └── public
│   │   ├── external_1 -> path/to/external/data
│       └── internal
```

The user should generate an `R/external_path.R` file and add those paths for each public or private external data. The following is an example of the `R/external_path.R` file.

```sh
DATA_PRIVATE_EXT_1 <- "~/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/private/data_1"
DATA_PUBLIC_EXT_1 <- "~/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/public/data_1"
DATA_PRIVATE_EXT_2 <- "~/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/private/data_2"
DATA_PUBLIC_EXT_2 <- "~/mqin_pm_no2_o3-adrd_hosp-medicare-causalgps/public/data_2"
```

There can be many connections to external folders. However, there is one location for internal data. You can put internal data inside your internal directories and have access to them directly. The `path_obj` includes all path information you will need to use in your code to ensure the correct location for the project. These paths will be automatically set internally. So using `path_obj,` instead of hardcoding the path, makes it easier for other members to checkout your code and analyze but store the results in their folder project without interfering with your processing. A sub-project can have one or many connections. 

- **Note**: As one can see, we do not keep the used external paths (`R/external_path.R`) in the committed code. The reason for this is to separate hard-coded paths from code as much as possible. If you wish to regenerate the data, you can see the content of `R/external_path.R` in the `sessioninfo` file inside the project's results folder. 

### Body of the research

This section is all codes that you would include if there was no convention. There are a couple of notes to consider:

- We do not use `library()` to import libraries. Instead of that, add the library name into the `DESCRIPTION` file and use the library name before their command. 
  - For example, `data.table::setDT(mydata)`
  - After any changes to `DESCRIPTION` run:
    - `devtools::document()`
    - `devtools::load_all()`

- Use `path_obj` to have access to all input and output files and folders. Available paths include:
   - `sp_dir`: An upper-level output directory for sub-project
   - `sp_output`: Output directory for sub-project (this includes data and figures)
   - `sp_cache`: Directory to hold cache values for the sub-project. The user should be able to remove the contents of this folder without any problem. The content of this folder is automatically generated and used. It is just a place to prevent computing the same internal results over and over again.
   - `sp_log`: Directory to host log files. 
   - `dir_data_private_ext_1`: Directory that points to your provided external private data directory. 
   - `private_int_ddir`: Directory that points to the internal private data directory. The actual data files should be in this directory.
   - `dir_data_private_ext_1`: Directory that points to your provided external public data directory. 
   - `public_int_ddir`: Directory that points to internal public data directory. The actual data files should be in this directory.

```r
# [Fourth section]
# Body of the research
#
# Put all your scientific analyses codes here.
#
#

# example of storing figure

pdf(file.path(path_obj$sp_output,"density_plot.pdf"), width = 8, height = 8)
plot(density(data$tmmx))
dev.off()
```
### Sub-project termination

This is an important step in documenting the environment. This includes recording the session info as well as the hash value of the code that is used. Later on, you will know which code version, specifically which commit, was used during the process. 

```r
# [Fifth section]
# Termination

finalize_subproject(path_obj = path_obj)
```
Anytime you run the code for each sub-project, a new session info file will be generated; the file name has the date and time of generation, so you can sort them chronologically. 


## Other Topics

- Any processing results will be in the `results` folder. Please note that the content of the results folder is never pushed to git. If you wish to share your results, you need to do it manually. 
- The `NEWS.md` is a good place to add your recent changes and any development with the project. 

