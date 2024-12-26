###=========================================================================#
### dalymod settings // inspired by knitr 
###=========================================================================#

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- settings ........... list object containing all settings functions
###---| settings_get ..... get settings by name
###---| settings_set ..... set settings
###-----| resolve ........ extract ... elements as list
###-----| merge .......... merge new settings into existing list
###-------| merge_list ... merge two lists with overlapping names

settings_get <-
  function(name) {
    dalysettings[name]
  }

settings_set <-
  function(...) {
    dalysettings <<- merge(resolve(...))
  }

resolve <- function(...) {
  dots <- list(...)
  if (length(dots) == 0) return()
  if (is.null(names(dots)) && length(dots) == 1 && is.list(dots[[1]]))
    if (length(dots <- dots[[1]]) == 0) return()
  return(dots)
}

merge <-
  function(values) {
    merge_list(dalysettings, values)
  }

merge_list <-
  function(x, y) {
    x[names(y)] <- y
    x
  }

settings <-
  list(
    get = settings_get,
    set = settings_set)
