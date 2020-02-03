# pipeline development coding 

# GOAL 1: extract function names and functions from a package

# It turns out Goal 1 is pretty easy. 
# You can use the function "getNamespaceExports" to get all the exported functions
# from a package, returned as a vector. 

package_name <- "x3ptools"
export_list <- getNamespaceExports(package_name) # from base R

# Then, you can grab any function from that vector and get the function itself
# by using the `lazyeval` package, specifically the `lazy_eval` function. 


function_chr <- paste0(package_name, "::", export_list[25])
function_text <- lazyeval::lazy_eval(function_chr) # from lazyeval package 

# GOAL 2: extract metadata about package including: 
# CRAN or GitHub
# CRAN version and GitHub most recent commit

package_information <- devtools::package_info(package_name) # from devtools package
package_information[package_name,]$source 
# for x3ptools, gives "Github (heike/x3ptools@d3ac438)", 
# so we can use this to get Github and most recent commit
# for ggplot2, gives "CRAN (R 3.6.0)", and if CRAN is true, 
package_information[package_name,]$loadedversion
# gives the version used 


# GOAL 3: Identify all the functions used in a file 

# It turns out there are lots of bread crumbs for this,
# and I think it should be relatively easy to pull out all the functions used. 
# The problem is identifying the package each function is called from.
filename <- "~/Desktop/IowaState/CSAFE/CSAFE-ISU/bullet-scan-variability/development/data_extract_by_barrel_round.R"
file_parsed <- getParseData(parse(filename)) # getParseData is from the utils package
functions_used <- file_parsed[file_parsed$token == "SYMBOL_FUNCTION_CALL",]$text

# interesting note about getParseData; 
# the following doesn't identify x3p_m_to_mum as a SYMBOL_FUNCTION_CALL, just : 
#scans <- scans %>% mutate(
#  x3p = x3p %>% purrr::map(.f = x3p_m_to_mum)
#)

# but this one does identify x3p_to_df as a SYMBOL_FUNCTION_CALL: 
#scans <- scans %>% mutate(ccdata = purrr::map(x3p, .f = function(x3p){
#  x3ptools::x3p_to_df(x3p)
#}))

# We should be able to take a list of the functions in a package that is used, 
# Check the script for each function, and subset on only the functions that are used, 
# And add that to the list. 

# So, this one still needs some work.  
# However, when we take the list of functions used, we can get the environment it is from.
# A major note, however, is that I'm unsure how this would work if there are multiple packages loaded which have the function.
# As far as I can tell, it will default to the most recent package that was loaded. 
sapply(functions_used, FUN = function(x) environmentName(environment(get(x))))
sapply(functions_used, FUN = function(x) environmentName(environment(match.fun(x)))) 
# These two approaches are identical with my small example. 

# Let's do a small example: 
# currently, "select" is in the "dplyr" environment.
library(dplyr)
environmentName(environment(select))
library(MASS)
environmentName(environment(select))
# Now it is in the "MASS" environment, because of how things are ordered in the environment. 
# One option is to use the conflicted package! 


# GOAL 4: Take a list of packages called in a script, generate list of functions, 
# and search for them in that script. 

# get rows where library is called:
library_calls <- readLines(filename)[str_detect(readLines(filename), "library") == T]
package_calls <- str_sub(library_calls, 9, -2)

package_meta <- data.frame(package_calls = package_calls, stringsAsFactors = F)
package_meta <- package_meta %>% mutate(package_exports = package_calls %>% purrr::map(.f = function(package){
  getNamespaceExports(package)
}))

filename <- "~/Desktop/IowaState/CSAFE/CSAFE-ISU/bullet-scan-variability/development/data_extract_by_barrel_round.R"
file_parsed <- getParseData(parse(filename))

package_meta <- package_meta %>% mutate(functions_called = package_exports %>% purrr::map(.f = function(package_exports){
  file_parsed <- file_parsed <- getParseData(parse(filename))
  sapply(package_exports, FUN = function(x) {sum(file_parsed$text == x)})
}))


# issue with "tidyverse" is that it's Namespace Exports aren't the functions in all the packages...



# Some general things I found when searching: 

# checkpoint package: checkpoint saves a snapshot of CRAN every day, 
# so when collaborating, you can say "checkpoint(date)" and it will use 
# everything as it was on CRAN at that date! Kind of cool. 

# conflicted package: throws an error when you call a function from multiple packages in your environment
# I think this may be a good thing to incorporate as a catch for what we are working on. 
# example from the repo: 
# filter(mtcars, cyl == 8)
#> Error: [conflicted] `filter` found in 2 packages.
#> Either pick the one you want with `::` 
#> * dplyr::filter
#> * stats::filter
#> Or declare a preference with `conflicted_prefer()`
#> * conflict_prefer("filter", "dplyr")
#> * conflict_prefer("filter", "stats")

# pipeR package: allows you to chain commands into a "pipeline"
#e.g. pipeline(1:10, mean, sum)

# dials package: tuning parameters and stuff
