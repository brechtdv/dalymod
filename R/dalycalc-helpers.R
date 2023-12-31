###=========================================================================#
### dalycalc helpers
###=========================================================================#

library(dplyr)

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- dalycalc .............. main wrapper
###-- dalycalc .............. main wrapper

split_age_string <-
  function(age) {
    age_full <-
      if (grepl("<", age)) {
        seq(0, as.numeric(gsub("<", "", age)), 5)
        
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
        sex_full
      }
    return(sex_full)
  }

split_agesex <-
  function(agesex_agg, country, pop) {
    agesex_agg_pop <-
      expand.grid(
        AGE = split_age_string(agesex_agg$AGE),
        SEX = split_sex_string(agesex_agg$SEX))
    agesex_agg_pop <-
      merge(agesex_agg_pop, subset(pop, ISO3 == country))
    agesex_agg_pop$W <-
      agesex_agg_pop$POP / sum(agesex_agg_pop$POP)
    agesex_agg_pop$INC_NR <-
      lapply(agesex_agg_pop$W, function(x) x * agesex_agg$INC_NR)
    
    agesex_agg_pop$W <- NULL
    agesex_agg_pop$ISO3 <- NULL
    
    apply(agesex_agg_pop, 1, as.list)
  }

split_agesex_all <-
  function(agesex_agg_all, country, pop) {
    agesex_all <- sapply(agesex_agg_all, split_agesex, country, pop)
    agesex_all <- unlist(agesex_all, recursive = FALSE)
    return(agesex_all)
  }

## main wrapper for DALY calculations
## .. requires processed 'dalymod' object as input
## .. 'year' is single numerical value
## .. 'pop' is dataframe with columns 'POP', 'ISO3', 'AGE', 'SEX'
## .. 'rle' is dataframe with columns 'AGE', 'SEX', 'RLE'
dalycalc <-
  function(.dalymod, year, pop, rle) {
    ## aggregate pop by country
    pop_country <- aggregate(POP ~ ISO3, pop, sum)
    
    ## calculate DALYs across YLD/YLL nodes
    ## .. .dalycalc is list of node>country>agesex
    nodes <- subset(.dalymod$dismod, CONTRIBUTION != "NA")$NODE
    .dalycalc <-
      lapply(.dalymod$nodes[nodes],
        dalycalc_node, year, pop_country, pop, rle)
    
    ## return dalycalc list
    return(.dalycalc)
  }

dalycalc_node <-
  function(node, year, pop_country, pop, rle) {
    ## expand countries if needed
    if (node$inc$COUNTRY == "ALL") {
      node_inc <- cbind(pop_country$ISO3, node$inc)
      node_inc$COUNTRY <- NULL
      names(node_inc)[names(node_inc) == "pop_country$ISO3"] <- "COUNTRY"
    }

    ## add year if needed
    if (all(node_inc$YEAR == "ALL")) {
      node_inc$YEAR <- year
    } else {
      node_inc <- subset(node_inc, YEAR == year)
    }
    
    ## merge incidence and population
    node_inc <- merge(node_inc, pop_country, by.x = "COUNTRY", by.y = "ISO3")
    
    ## multiply inc with pop
    node_inc$INC_NR <-
      apply(node_inc, 1,
            function(x) x[["POP"]] * x[["SAMPLES"]], simplify = FALSE)
    
    ## setup age split
    age_split <- node$age$samp
    age_split$YEAR <- NULL
    if (age_split$COUNTRY == "ALL") {
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
    if (sex_split$COUNTRY == "ALL") {
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
      if (node$dur$samp$COUNTRY == "ALL") {
        node_dur <- cbind(pop_country$ISO3, node$dur$samp)
        node_dur$COUNTRY <- NULL
        names(node_dur)[names(node_dur) == "pop_country$ISO3"] <- "COUNTRY"
      } else {
        node_dur <- node$dur$samp
      }
      if (node$dsw$samp$COUNTRY == "ALL") {
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
      
      for (i in seq_along(dalycalc_all)) { # country
        for (j in seq_along(dalycalc_all[[i]])) { # age-sex
          dalycalc_all[[i]][[j]]$YLD_NR <-
            dalycalc_all[[i]][[j]]$INC_NR *
              node_dur$SAMPLES[[i]] *
              node_dsw$SAMPLES[[i]]
        }
      }
    }

    ## calculate YLL
    if (node$set$contribution == "YLL") {
      for (i in seq_along(dalycalc_all)) { # country
        for (j in seq_along(dalycalc_all[[i]])){ # age-sex
          dalycalc_all[[i]][[j]]$YLL_NR <-
            dalycalc_all[[i]][[j]]$INC_NR * 
              get_rle(
                dalycalc_all[[i]][[j]]$AGE,
                dalycalc_all[[i]][[j]]$SEX,
                rle)
        }
      }
    }
    
    ## return dalycalc node
    return(dalycalc_all)
}

get_rle <-
  function(age, sex, rle) {
    subset(rle, AGE == age & SEX == sex)$RLE
  }

##


## import
source("R/dalymod-helpers.R")
file <-
  "C:/Users/BrDe394/Dropbox/#FERG2/#CTF/daly-calculation/ferg2-daly-anthrax.xlsx"
.dalymod <- dalymod(file, 5)

## set year
year <- 2010

## aggregate population by age groups
pop <- subset(FERG2:::pop, YEAR == year)
pop$AGE2 <- cut(pop$AGE, c(seq(0, 85, 5), Inf), right = FALSE)
pop2 <- aggregate(POP ~ AGE2 + SEX + ISO3, pop, sum)
names(pop2)[names(pop2) == "AGE2"] <- "AGE"
str(pop2)
pop <- pop2

## residual life expectancy
# rle <- FERG2:::RLE
rle <- cbind(pop2[1:36, 1:2], RLE = 1:36)


## this creates the detailed calculations
## .. all nodes, all countries, all age/sex groups
.dalycalc <- dalycalc(.dalymod, year, pop, rle)
saveRDS(.dalycalc, file = "dalycalc-anthrax.rds")

## this aggregates estimates across nodes
## .. all countries, all age/sex groups
.dalycalc_country <- dalycalc_aggregate_nodes(.dalycalc)
str(.dalycalc_country)
.dalycalc_country$BEL
.dalycalc_country_age$BEL
.dalycalc_country$BEL[[1]]

agesex <-
  list(
    "child" =
      data.frame(
        AGE = c("[0,5)", "[0,5)"),
        SEX = c("Male", "Female")),
    "adult" =
      data.frame(
        AGE = c("[5,10)", "[5,10)", "[10,15)", "[10,15)"),
        SEX = c("Male", "Female", "Male", "Female")))

## this aggregates estimates across age groups
## .. all countries, broad age/sex groups
.dalycalc_country_age <-
  dalycalc_aggregate_agesex(.dalycalc_country, agesex)
str(.dalycalc_country_age)
.dalycalc_country_age[["BEL"]][[1]]
.dalycalc_country_age[["BEL"]][[2]]

## TO DO

## .. export INC contrib
## .. calculate aggregates across countries
## .. export samples > dalycalc_samples

list_sum <-
  function(x) {
    rowSums(do.call("cbind", x))
  }

dalycalc_aggregate_nodes <-
  function(.dalycalc) {
    dalycalc_agg <- vector("list", length(.dalycalc[[1]]))
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
        dalycalc_agg[[i]][[j]]$YLD_NR <-
          list_sum(lapply(.dalycalc, function(x) x[[i]][[j]]$YLD_NR))
        dalycalc_agg[[i]][[j]]$YLL_NR <-
          list_sum(lapply(.dalycalc, function(x) x[[i]][[j]]$YLL_NR))
        dalycalc_agg[[i]][[j]]$DALY_NR <-
          dalycalc_agg[[i]][[j]]$YLD_NR + dalycalc_agg[[i]][[j]]$YLL_NR
      }
    }
    
    return(dalycalc_agg)
  }

dalycalc_aggregate_agesex <-
  function(.dalycalc, agesex) {
    agesex_names <- names(agesex)
    age_names <- sapply(.dalycalc[[1]], function(x) x$AGE)
    sex_names <- sapply(.dalycalc[[1]], function(x) x$SEX)
    agesex_all <- data.frame(AGE = age_names, SEX = sex_names)
    
    dalycalc_agg <- vector("list", length(.dalycalc))
    names(dalycalc_agg) <- names(.dalycalc)
      
    for (i in seq_along(dalycalc_agg)) {
      dalycalc_agg[[i]] <- vector("list", length(agesex_names))
      names(dalycalc_agg[[i]]) <- agesex_names
      for (j in seq_along(agesex_names)) { # age-sex
        agesex_id <-
          which(duplicated(rbind(agesex_all, agesex[[j]]), fromLast = TRUE))
        dalycalc_agg[[i]][[j]]$POP <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$POP))
        dalycalc_agg[[i]][[j]]$INC_NR <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$INC_NR))
        dalycalc_agg[[i]][[j]]$YLD_NR <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$YLD_NR))
        dalycalc_agg[[i]][[j]]$YLL_NR <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$YLL_NR))
        dalycalc_agg[[i]][[j]]$DALY_NR <-
          list_sum(lapply(.dalycalc[[i]][agesex_id], function(x) x$DALY_NR))
      }
    }

    return(dalycalc_agg)
  }

.dalycalc_age <- .dalycalc_country_age

## get samples in long format
dalycalc_samples <-
  function() {
    
  }

?geom_density

## this should work on .dalycalc_country and .dalycalc_country_age
dalycalc_summary <-
  function(.dalycalc_agg, pars = c("INC_NR", "YLD_NR", "YLL_NR'")) {
    names(pars) <- pars
    out_lst <- lapply(pars, dalycalc_summary_par, .dalycalc_agg = .dalycalc_agg)
    dplyr::bind_rows(out_lst, .id = "MEASURE")
    return(out_df)
  }
    
dalycalc_summary_par <-
  function(.dalycalc_agg, par) {
    out_lst <- lapply(.dalycalc_agg, dalycalc_summary_par_age, par)
    dplyr::bind_rows(out_lst, .id = "COUNTRY")
    return(out_df)
  }

dalycalc_summary_par_age <-
  function(.dalycalc_agg_age, par) {
    out <-
    data.frame(
      AGE = sapply(.dalycalc_agg_age[[1]], function(x) x$AGE),
      SEX = sapply(.dalycalc_agg_age[[1]], function(x) x$SEX),
      t(sapply(.dalycalc_agg_age[[1]], function(x) bd::mean_ci(x[[par]]))))
    names(out) <- c("AGE", "SEX", "VAL_MEAN", "VAL_LWR", "VAL_UPR")
    return(out)
  }

dalycalc_summary(.dalycalc_country_age)
dalycalc_summary(.dalycalc_country_age)
## does not work because 'AGE' and 'SEX' not found


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


##
##

file <-
  "C:/Users/BrDe394/Dropbox/#FERG2/#CTF/daly-calculation/ferg2-daly-anthrax.xlsx"
dalymods <- dalymod(file, 5)
object.size(dalymods)

