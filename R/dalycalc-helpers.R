###=========================================================================#
### dalycalc helpers
###=========================================================================#

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- split_agesex_all ............... main wrapper
###---| split_agesex ................. main wrapper
###-----| split_age_string ........... main wrapper
###-----| split_sex_string ........... main wrapper
###-- get_rle ........................ main wrapper
###-- list_sum ....................... main wrapper
###-- dalycalc ....................... main wrapper
###---| dalycalc_node ................ main wrapper
###-- dalycalc_aggregate_nodes ....... main wrapper
###-- dalycalc_aggregate_agesex ...... main wrapper
###-- dalycalc_aggregate_country ..... main wrapper
###-- dalycalc_summary ............... main wrapper
###---| dalycalc_summary_par ......... main wrapper
###-----| dalycalc_summary_par_age ... main wrapper
###-- dalycalc_add ................... main wrapper
###-- dalycalc_mult .................. main wrapper

##--------------------------------------------------------------------------#
## generic helpers ---------------------------------------------------------#

split_agesex_all <-
  function(agesex_agg_all, country, pop) {
    agesex_all <- lapply(agesex_agg_all, split_agesex, country, pop)
    agesex_all <- unlist(agesex_all, recursive = FALSE)
    return(agesex_all)
  }

split_agesex <-
  function(agesex_agg, country, pop) {
    agesex_agg_pop <-
      expand.grid(
        AGE = split_age_string(agesex_agg$AGE),
        SEX = split_sex_string(agesex_agg$SEX))
    agesex_agg_pop <-
      base::merge(agesex_agg_pop, subset(pop, ISO3 == country))
    agesex_agg_pop$W <-
      agesex_agg_pop$POP / sum(agesex_agg_pop$POP)
    # replace NaN weight by zero
    agesex_agg_pop$W[is.nan(agesex_agg_pop$W)] <- 0
    agesex_agg_pop$INC_NR <-
      lapply(agesex_agg_pop$W, function(x) x * agesex_agg$INC_NR)
    
    agesex_agg_pop$W <- NULL
    agesex_agg_pop$ISO3 <- NULL
    
    apply(agesex_agg_pop, 1, as.list)
  }

split_age_string <-
  function(age) {
    age_full <-
      if (age == "<1") {
        c(0, 1)
        
      } else if (age == "1-4") {
        c(1, 5)
        
      } else if (age == "1+") {
        c(1, seq(5, 85, 5), Inf)
        
      } else if (age == "0+") {
        c(0, 1, seq(5, 85, 5), Inf)
        
      } else if (grepl("<", age)) {
        c(0, 1, seq(5, as.numeric(gsub("<", "", age)), 5))
        
      } else if (grepl("\\+", age)) {
        c(seq(as.numeric(gsub("\\+", "", age)), 85, 5), Inf)
        
      } else {
        seq(as.numeric(gsub("-.*", "", age)),
            as.numeric(gsub(".*-", "", age)) + 1,
            5)
      }
    age_full <- levels(cut(1:100, age_full, right = FALSE))
    return(age_full)
  }

split_sex_string <-
  function(sex) {
    sex_full <-
      if (sex == "Both sexes") {
        c("Male", "Female")
        
      } else {
        sex
      }
    return(sex_full)
  }

get_rle <-
  function(age, sex, rle) {
    subset(rle, AGE == age & SEX == sex)$RLE
  }

list_sum <-
  function(x) {
    # compile results as matrix
    y <- do.call("cbind", x)
    if (is.null(y)) { # if empty, return 0
      0
    } else { # if not empty, calculate row sums
      rowSums(y)
    }
  }

## main wrapper for DALY calculations
## .. requires processed 'dalymod' object as input
## .. 'year' is single numerical value
## .. 'pop' is dataframe with columns 'POP', 'ISO3', 'AGE', 'SEX'
## .. 'rle' is dataframe with columns 'AGE', 'SEX', 'RLE'
## .. 'lle' is dataframe with columns 'AGE', 'SEX', 'LLE'
dalycalc <-
  function(.dalymod, year, pop, rle, lle = NULL, verbose = TRUE) {
    ## set cli 'cli.default_handler' option
    if (verbose) {
      options(cli.default_handler = NULL)
    } else {
      options(cli.default_handler = function(...) { })
    }
    
    ## aggregate pop by country
    # pop_country <- aggregate(POP ~ ISO3, pop, sum)
    
    ## calculate DALYs across YLD/YLL nodes
    ## .. .dalycalc is list of node>country>agesex
    nodes <-
      subset(.dalymod$dismod, CONTRIBUTION != "NA" | INCIDENCE == "YES")$NODE
    cli_progress_step("Calculating DALYs across nodes", spinner = TRUE)
    # .dalycalc <-
    #   lapply(.dalymod$nodes[nodes],
    #          dalycalc_node, year, pop_country, pop, rle)
    .dalycalc <- lapply(.dalymod$nodes[nodes], dalycalc_node, year, pop, rle)
    
    ## add S3 class
    class(.dalycalc) <- "dalycalc"
      
    ## return dalycalc list
    return(.dalycalc)
  }

dalycalc_node <-
  # function(node, year, pop_country, pop, rle) {
  function(node, year, pop, rle) {
    ## aggregate population by country
    ## .. this calculates correct pop denominator for age-specific nodes
    if (grepl("#", node$set$name)) {
      node_age <- gsub(".*#", "", node$set$name)
      node_age <- split_age_string(node_age)
      pop_country <- aggregate(POP ~ ISO3, subset(pop, AGE %in% node_age), sum)
      
    } else {
      pop_country <- aggregate(POP ~ ISO3, pop, sum)
    }
    
    ## extract incidence
    node_inc <- node$inc
    
    ## add or subset year if needed
    if (all(node_inc$YEAR == "ALL")) {
      node_inc$YEAR <- year
    } else {
      node_inc <- subset(node_inc, YEAR == year)
    }
    
    ## expand countries if needed
    if (all(node_inc$COUNTRY == "ALL")) {
      node_inc <- cbind(pop_country$ISO3, node$inc)
      node_inc$COUNTRY <- NULL
      names(node_inc)[names(node_inc) == "pop_country$ISO3"] <- "COUNTRY"
    }
    
    ## merge incidence and population
    node_inc <-
      base::merge(node_inc, pop_country, by.x = "COUNTRY", by.y = "ISO3")
    
    ## multiply inc with pop
    node_inc$INC_NR <-
      apply(node_inc, 1,
            function(x) x[["POP"]] * x[["SAMPLES"]], simplify = FALSE)
    
    ## setup age split
    age_split <- node$age$samp
    age_split$YEAR <- NULL
    if (all(age_split$COUNTRY == "ALL")) {
      age_split <- cbind(pop_country$ISO3, age_split)
      age_split$COUNTRY <- NULL
      names(age_split)[names(age_split) == "pop_country$ISO3"] <- "COUNTRY"
    }
    rownames(age_split) <- age_split$COUNTRY
    age_split$COUNTRY <- NULL
    age_names <- names(age_split)
    
    ## setup sex split
    sex_split <- node$sex$samp
    sex_split$YEAR <- NULL
    if (all(sex_split$COUNTRY == "ALL")) {
      sex_split <- cbind(pop_country$ISO3, sex_split)
      sex_split$COUNTRY <- NULL
      names(sex_split)[names(sex_split) == "pop_country$ISO3"] <- "COUNTRY"
    }
    rownames(sex_split) <- sex_split$COUNTRY
    sex_split$COUNTRY <- NULL
    sex_names <- names(sex_split)
    
    ## stratify cases by age ~ dist
    dalycalc_agg <- vector("list", nrow(node_inc))
    names(dalycalc_agg) <- node_inc$COUNTRY

    # loop over countries
    for (i in seq_along(dalycalc_agg)) {
      id <- 0
      dalycalc_agg[[i]] <-
        vector("list", length(age_names) * length(sex_names))
      
      # loop over age groups
      for (j in seq_along(age_names)) {
      
        # loop over sexes
        for (k in seq_along(sex_names)) {
          id <- id + 1
          dalycalc_agg[[i]][[id]]$AGE <- age_names[j]
          dalycalc_agg[[i]][[id]]$SEX <- sex_names[k]
          dalycalc_agg[[i]][[id]]$INC_NR <-
            node_inc$INC_NR[[i]] *
              age_split[[age_names[j]]][[i]] *
              sex_split[[sex_names[k]]][[i]]
        }
      }
    }

    ## stratify cases by age ~ full
    dalycalc_all <-
    lapply(
      names(dalycalc_agg),
      function(x) split_agesex_all(dalycalc_agg[[x]], x, pop))
    names(dalycalc_all) <- names(dalycalc_agg)

    ## calculate YLD
    if (node$set$contribution == "YLD") {
      
      ## expand countries if needed
      if (all(node$dur$samp$COUNTRY == "ALL")) {
        node_dur <- cbind(pop_country$ISO3, node$dur$samp)
        node_dur$COUNTRY <- NULL
        names(node_dur)[names(node_dur) == "pop_country$ISO3"] <- "COUNTRY"
      } else {
        node_dur <- node$dur$samp
      }
      if (all(node$dsw$samp$COUNTRY == "ALL")) {
        node_dsw <- cbind(pop_country$ISO3, node$dsw$samp)
        node_dsw$COUNTRY <- NULL
        names(node_dsw)[names(node_dsw) == "pop_country$ISO3"] <- "COUNTRY"
      } else {
        node_dsw <- node$dsw$samp
      }
      
      ## add year if needed
      if (all(node_dur$YEAR == "ALL")) {
        node_dur$YEAR <- year
      } else {
        node_dur <- subset(node_dur, YEAR == year)
      }
      if (all(node_dsw$YEAR == "ALL")) {
        node_dsw$YEAR <- year
      } else {
        node_dsw <- subset(node_dsw, YEAR == year)
      }
      
      # if lifelong duration, get duration from 'lle' dataframe
      # otherwise, use duration from dalymod
      # TO DO > make LLE sex specific!
      
      if ("LIFELONG" %in% names(node_dur)) {
        for (i in seq_along(dalycalc_all)) { # country
          for (j in seq_along(dalycalc_all[[i]])) { # age-sex
            dalycalc_all[[i]][[j]]$YLD_NR <-
              dalycalc_all[[i]][[j]]$INC_NR *
              subset(
                lle,
                ISO3 == names(dalycalc_all)[i] &
                  AGE == dalycalc_all[[i]][[j]]$AGE &
                  #SEX == dalycalc_all[[i]][[j]]$SEX &
                  YEAR == dalycalc_all[[i]][[j]]$YEAR)$LLE *
              node_dsw$SAMPLES[[i]]
          }
        }
        
      } else {
        for (i in seq_along(dalycalc_all)) { # country
          for (j in seq_along(dalycalc_all[[i]])) { # age-sex
            dalycalc_all[[i]][[j]]$YLD_NR <-
              dalycalc_all[[i]][[j]]$INC_NR *
              node_dur$SAMPLES[[i]] *
              node_dsw$SAMPLES[[i]]
          }
        }
      }
      
      ## set INC_NR to zero if no INC contribution
      if (node$set$incidence == "NO") {
        for (i in seq_along(dalycalc_all)) { # country
          for (j in seq_along(dalycalc_all[[i]])){ # age-sex
            dalycalc_all[[i]][[j]]$INC_NR <- NULL
          }
        }
      }

    }
    
    ## calculate YLL & rename INC_NR to MRT_NR
    if (node$set$contribution == "YLL") {
      for (i in seq_along(dalycalc_all)) { # country
        for (j in seq_along(dalycalc_all[[i]])){ # age-sex
          dalycalc_all[[i]][[j]]$YLL_NR <-
            dalycalc_all[[i]][[j]]$INC_NR * 
              get_rle(
                dalycalc_all[[i]][[j]]$AGE,
                dalycalc_all[[i]][[j]]$SEX,
                rle)
          dalycalc_all[[i]][[j]]$MRT_NR <-
            dalycalc_all[[i]][[j]]$INC_NR
          dalycalc_all[[i]][[j]]$INC_NR <- NULL
        }
      }
    }
    
    ## impose consistent order
    order_to <-
      paste(
        rep(split_age_string("0+"), each = 2),
        rep(c("Male", "Female"), 19))
    order_from <-
      paste(
        sapply(dalycalc_all[[1]], function(x) x$AGE),
        sapply(dalycalc_all[[1]], function(x) x$SEX))
    for (i in seq_along(dalycalc_all)) { # country
      dalycalc_all[[i]] <- dalycalc_all[[i]][match(order_from, order_to)]
    }
    
    ## return dalycalc node
    return(dalycalc_all)
}

## TO DO

## .. export INC contrib
## .. export samples > dalycalc_samples

dalycalc_aggregate_nodes <-
  function(.dalycalc) {
    ## check class
    if (!("dalycalc" %in% class(.dalycalc)))
      stop("Input must be of class ", sQuote("dalycalc"))
    
    ## aggregate nodes
    dalycalc_agg <- vector("list", length(.dalycalc[[1]])) # country
    names(dalycalc_agg) <- names(.dalycalc[[1]])
    age_names <- sapply(.dalycalc[[1]][[1]], function(x) x$AGE)
    sex_names <- sapply(.dalycalc[[1]][[1]], function(x) x$SEX)
    
    for (i in seq_along(dalycalc_agg)) { # country
      dalycalc_agg[[i]] <- vector("list", length(age_names))
      pop <- sapply(.dalycalc[[1]][[i]], function(x) x$POP)
      for (j in seq_along(age_names)) { # age-sex
        dalycalc_agg[[i]][[j]]$AGE <- age_names[j]
        dalycalc_agg[[i]][[j]]$SEX <- sex_names[j]
        dalycalc_agg[[i]][[j]]$POP <- pop[j]
        dalycalc_agg[[i]][[j]]$INC_NR <-
          list_sum(lapply(.dalycalc, function(x) x[[i]][[j]]$INC_NR))
        dalycalc_agg[[i]][[j]]$MRT_NR <-
          list_sum(lapply(.dalycalc, function(x) x[[i]][[j]]$MRT_NR))
        dalycalc_agg[[i]][[j]]$YLD_NR <-
          list_sum(lapply(.dalycalc, function(x) x[[i]][[j]]$YLD_NR))
        dalycalc_agg[[i]][[j]]$YLL_NR <-
          list_sum(lapply(.dalycalc, function(x) x[[i]][[j]]$YLL_NR))
        dalycalc_agg[[i]][[j]]$DALY_NR <-
          dalycalc_agg[[i]][[j]]$YLD_NR + dalycalc_agg[[i]][[j]]$YLL_NR
      }
    }
    
    ## add S3 class
    class(dalycalc_agg) <- "dalycalc_agg"
    
    ## return output
    return(dalycalc_agg)
  }

dalycalc_aggregate_agesex <-
  function(.dalycalc, age = NULL, sex = NULL) {
    ## check class
    if (!("dalycalc_agg" %in% class(.dalycalc)))
      stop("Input must be of class ", sQuote("dalycalc_agg"))
    
    ## get age and sex names from old definition
    age_old <- sapply(.dalycalc[[1]], function(x) x$AGE)
    sex_old <- sapply(.dalycalc[[1]], function(x) x$SEX)
    agesex_old <- data.frame(AGE = age_old, SEX = sex_old)
    
    ## if age or sex not specified, reuse old definition
    if (is.null(age)) {
      age <- as.list(unique(age_old))
      names(age) <- unique(age_old)
    }
    
    if (is.null(sex)) {
      sex <- as.list(unique(sex_old))
      names(sex) <- unique(sex_old)
    }

    ## compile new age and sex definitions 
    age_df <-
      data.frame(
        AGE = unlist(age),
        AGE_NEW = rep(names(age), times = sapply(age, length)))
    sex_df <-
      data.frame(
        SEX = unlist(sex),
        SEX_NEW = rep(names(sex), times = sapply(sex, length)))
    agesex_new <-
      expand.grid(AGE = names(age), SEX = names(sex))
    
    ## merge old and new definitions
    agesex_oldnew <-
      base::merge(base::merge(agesex_old, age_df), sex_df)
    
    ## restore order of 'agesex_old'
    agesex_oldnew <- base::merge(agesex_old, agesex_oldnew, sort = FALSE)
    
    ## prepare object to hold results
    dalycalc_agg <- vector("list", length(.dalycalc))
    names(dalycalc_agg) <- names(.dalycalc)
      
    ## add up numbers across country and new age-sex group
    for (i in seq_along(dalycalc_agg)) {
      dalycalc_agg[[i]] <- vector("list", nrow(agesex_new))
      for (j in seq(nrow(agesex_new))) { # age-sex
        # identify rows in old age-sex definition
        agesex_id <-
          with(agesex_oldnew,
               which(AGE_NEW == agesex_new[j, "AGE"] &
                     SEX_NEW == agesex_new[j, "SEX"]))
        
        # use new age-sex definition
        dalycalc_agg[[i]][[j]]$AGE <- agesex_new[j, "AGE"]
        dalycalc_agg[[i]][[j]]$SEX <- agesex_new[j, "SEX"]
        
        # add up numbers
        dalycalc_agg[[i]][[j]]$POP <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$POP))
        dalycalc_agg[[i]][[j]]$INC_NR <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$INC_NR))
        dalycalc_agg[[i]][[j]]$MRT_NR <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$MRT_NR))
        dalycalc_agg[[i]][[j]]$YLD_NR <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$YLD_NR))
        dalycalc_agg[[i]][[j]]$YLL_NR <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$YLL_NR))
        dalycalc_agg[[i]][[j]]$DALY_NR <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$DALY_NR))
      }
    }
    
    ## add S3 class
    class(dalycalc_agg) <- "dalycalc_agg"
    
    ## return output
    return(dalycalc_agg)
  }

## aggregate results into regional groupings
dalycalc_aggregate_country <-
  function(.dalycalc, country) {
    ## check class
    if (!("dalycalc_agg" %in% class(.dalycalc)))
      stop("Input must be of class ", sQuote("dalycalc_agg"))
    
    ## prepare object to hold results
    dalycalc_agg <- vector("list", length(country))
    names(dalycalc_agg) <- names(country)
    
    ## add up numbers across region
    for (i in seq_along(dalycalc_agg)) { # new regions
      dalycalc_agg[[i]] <- vector("list", length(.dalycalc[[1]]))
      
      for (j in seq_along(.dalycalc[[1]])) { # age-sex
        # reapply age-sex definition
        dalycalc_agg[[i]][[j]]$AGE <- .dalycalc[[1]][[j]]$AGE
        dalycalc_agg[[i]][[j]]$SEX <- .dalycalc[[1]][[j]]$SEX
        
        # add up numbers
        dalycalc_agg[[i]][[j]]$POP <-
          list_sum(lapply(.dalycalc[country[[i]]], function(x) x[[j]]$POP))
        dalycalc_agg[[i]][[j]]$INC_NR <-
          list_sum(lapply(.dalycalc[country[[i]]], function(x) x[[j]]$INC_NR))
        dalycalc_agg[[i]][[j]]$MRT_NR <-
          list_sum(lapply(.dalycalc[country[[i]]], function(x) x[[j]]$MRT_NR))
        dalycalc_agg[[i]][[j]]$YLD_NR <-
          list_sum(lapply(.dalycalc[country[[i]]], function(x) x[[j]]$YLD_NR))
        dalycalc_agg[[i]][[j]]$YLL_NR <-
          list_sum(lapply(.dalycalc[country[[i]]], function(x) x[[j]]$YLL_NR))
        dalycalc_agg[[i]][[j]]$DALY_NR <-
          list_sum(lapply(.dalycalc[country[[i]]], function(x) x[[j]]$DALY_NR))
      }
    }
    
    ## add S3 class
    class(dalycalc_agg) <- "dalycalc_agg"
    
    ## return output
    return(dalycalc_agg)
  }

## get samples in long format
dalycalc_samples <-
  function() {
    
  }

## this should work on .dalycalc_country and .dalycalc_country_age
dalycalc_summary <-
  function(
    .dalycalc_agg,
    pars = c("INC_NR", "MRT_NR", "YLD_NR", "YLL_NR", "DALY_NR",
             "MRT_INC", "YLD_INC", "YLL_INC", "DALY_INC")) {
    ## check class
    if (!("dalycalc_agg" %in% class(.dalycalc_agg)))
      stop("Input must be of class ", sQuote("dalycalc_agg"))
    
    ## calculate summaries
    names(pars) <- pars
    out_lst <- lapply(pars, dalycalc_summary_par, .dalycalc_agg = .dalycalc_agg)
    out_df <- dplyr::bind_rows(out_lst, .id = "MEASURE")
    return(out_df)
  }

dalycalc_summary_par <-
  function(.dalycalc_agg, par) {
    out_lst <- lapply(.dalycalc_agg, dalycalc_summary_par_age, par)
    out_df <- dplyr::bind_rows(out_lst, .id = "COUNTRY")
    return(out_df)
  }

dalycalc_summary_par_age <-
  function(.dalycalc_agg_age, par) {
    # define numerator/denominator
    par1 <- paste0(gsub("_.*", "", par), "_NR")
    par2 <- paste0(gsub(".*_", "", par), "_NR")
    
    # add dummy denominator
    .dalycalc_agg_age <-
      lapply(.dalycalc_agg_age, 
             function(x) {
               x$NR_NR <- 1
               return(x)
             })
    
    # compile summaries
    out <-
      data.frame(
        AGE = sapply(.dalycalc_agg_age, function(x) x$AGE),
        SEX = sapply(.dalycalc_agg_age, function(x) x$SEX),
        POP = sapply(.dalycalc_agg_age, function(x) x$POP),
        VAL_MEAN =
          sapply(.dalycalc_agg_age,
                 function(x) mean(x[[par1]] / x[[par2]])),
        VAL_MEDIAN =
          sapply(.dalycalc_agg_age,
                 function(x) median(x[[par1]] / x[[par2]])),
        VAL_LWR =
          sapply(.dalycalc_agg_age,
                 function(x) quantile(x[[par1]] / x[[par2]], 0.025, na.rm = TRUE)),
        VAL_UPR =
          sapply(.dalycalc_agg_age,
                 function(x) quantile(x[[par1]] / x[[par2]], 0.975, na.rm = TRUE)))
    rownames(out) <- NULL
    return(out)
  }

dalycalc_add <-
  function(...) {
    x <- list(...)
    out <- x[[1]]
    out_names <- names(out)
    
    for (i in seq_along(x)[-1]) {  # objects
      
      if (!is.null(x[[i]])) {
        
        # check consistency in country names
        if (!all(out_names %in% names(x[[i]])))
          stop("Identified mismatch in country names.")
        if (!all(names(x[[i]]) %in% out_names))
          stop("Identified mismatch in country names.")
        
        for (j in out_names) {  # countries
          
          # check consistency in settings
          if (!all(sapply(out[[j]], function(x) x$POP) ==
                   sapply(x[[i]][[j]], function(x) x$POP)))
            stop("Identified mismatch in POP in country ", j)
          
          if (!all(sapply(out[[j]], function(x) x$AGE) ==
                   sapply(x[[i]][[j]], function(x) x$AGE)))
            stop("Identified mismatch in AGE in country ", j)
          
          if (!all(sapply(out[[j]], function(x) x$SEX) ==
                   sapply(x[[i]][[j]], function(x) x$SEX)))
            stop("Identified mismatch in SEX in country", j)
          
          for (k in seq_along(x[[i]][[j]])) {  # age-sex
            out[[j]][[k]]$INC_NR <- out[[j]][[k]]$INC_NR + x[[i]][[j]][[k]]$INC_NR
            out[[j]][[k]]$MRT_NR <- out[[j]][[k]]$MRT_NR + x[[i]][[j]][[k]]$MRT_NR
            out[[j]][[k]]$YLD_NR <- out[[j]][[k]]$YLD_NR + x[[i]][[j]][[k]]$YLD_NR
            out[[j]][[k]]$YLL_NR <- out[[j]][[k]]$YLL_NR + x[[i]][[j]][[k]]$YLL_NR
            out[[j]][[k]]$DALY_NR <- out[[j]][[k]]$DALY_NR + x[[i]][[j]][[k]]$DALY_NR
          }
        } 
      }
    }
    return(out) 
  }

dalycalc_mult <-
  function(.dalycalc, y) {
    ## check length match
    if (max(sapply(.dalycalc[[1]][[1]], length)) != length(y))
      warning("Lengths of objects do not match.")
    
    ## multiply nodes
    for (i in seq_along(.dalycalc)) {  # countries
      for (j in seq_along(.dalycalc[[i]])) {  # age-sex
        .dalycalc[[i]][[j]]$INC_NR <- .dalycalc[[i]][[j]]$INC_NR * y
        .dalycalc[[i]][[j]]$MRT_NR <- .dalycalc[[i]][[j]]$MRT_NR * y
        .dalycalc[[i]][[j]]$YLD_NR <- .dalycalc[[i]][[j]]$YLD_NR * y
        .dalycalc[[i]][[j]]$YLL_NR <- .dalycalc[[i]][[j]]$YLL_NR * y
        .dalycalc[[i]][[j]]$DALY_NR <- .dalycalc[[i]][[j]]$DALY_NR * y
      }
    }
    return(.dalycalc) 
  }

# dalycalc > per outcome, 5y age group
# COUNTRY/YEAR/AGE/SEX/POP/INC_RT/DUR/DSW/YLD_RT/INC_NR/YLD_NR
# COUNTRY/YEAR/AGE/SEX/POP/INC_RT/RLE/YLL_RT/INC_NR/YLL_NR
# 
# dalysum > grouped outcomes, 5vs5+
# COUNTRY/YEAR/AGE/SEX/POP/INC_RT/YLD_RT/INC_NR/YLD_NR
# COUNTRY/YEAR/AGE/SEX/POP/INC_RT/YLL_RT/INC_NR/YLL_NR
# 
# mean/UI
# COUNTRY/YEAR/AGE/SEX/METRIC/MEASURE/VAL_MEAN/VAL_LWR/VAL_UPR
# + regional groupings

