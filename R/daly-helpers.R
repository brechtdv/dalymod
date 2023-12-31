###=========================================================================#
### DALY calculation helpers
###=========================================================================#

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- read_excel_base ....... import excel file and return as 'data.frame'
###-- import_dalymod ........ import 'dalymod' excel and transform to list
###---| import_node ......... import 'dalymod' node and transform to list
###-- pre_sample_dalymod .... pre-sampler for dalymod object
###---| pre_sample_node ..... pre-sampler for one node
###-----| pre_sample_input .. pre-sampler for VAL,DUR,DSW
###-----| pre_sample_age .... pre-sampler for AGE type
###-----| pre_sample_sex .... pre-sampler for SEX type
###---------| sim_fixed ..... actual sampler, fixed
###---------| sim_beta ...... actual sampler, beta
###---------| sim_gamma ..... actual sampler, gamma
###---------| sim_unif ...... actual sampler, uniform
###---------| sim_pert ...... actual sampler, beta-pert
###---------| sim_mean_ci ... actual sampler, mean & 95%CI
###---------| optim_gamma ... find gamma pars through optimization
###---------| optim_beta .... find beta pars through optimization
###-- normalize_splits ...... normalize samples for splits
###-- multiply_dalymod
###---| multiply_nodes
###-----| get_tree
###-- dalymod
###-- dalycalc

library(readxl)
library(prevalence) #betaPERT

##--------------------------------------------------------------------------#
## generic helpers ---------------------------------------------------------#

read_excel_base <-
  function(...) {
    xl <- suppressMessages(read_excel(...))
    class(xl) <- "data.frame"
    return(xl)
  }

##--------------------------------------------------------------------------#
## import 'dalymod' excel and transform to list ----------------------------#

import_dalymod <-
  function(file) {
    ## check 'file' input
    if (!file.exists(file))
      stop("File ", sQuote(file), " does not exist.")
    if (!(tools::file_ext(file) %in% c("xls", "xlsx")))
      stop("File ", sQuote(file), " is not an Excel file.")
    
    ## read Excel sheets
    xl_sheets <- excel_sheets(file)
    if (!("DM" %in% xl_sheets))
      stop("Sheet 'DM' not found.")
    
    ## read Excel disease model
    xl_dismod <- read_excel_base(file, "DM")
    if (!all(xl_dismod$NODE %in% xl_sheets))
      stop("Not all disease model nodes are found as Excel sheets.")
    
    ## read Excel nodes
    nodes <- xl_dismod$NODE
    xl_nodes <- vector("list", length(nodes))
    names(xl_nodes) <- tolower(nodes)
    for (i in seq_along(nodes)) {
      xl_nodes[[i]] <- import_node(nodes[i])
    }
    
    ## return import
    return(list(dismod = xl_dismod, nodes = xl_nodes))
  }

import_node <-
  function(node) {
    ## read Excel node
    xl_node <- read_excel_base(file, node, col_names = FALSE)
    
    ## define row numbers
    row_type <- which(xl_node[[1]] == "TYPE")
    row_val <- which(xl_node[[1]] == "VALUE")
    row_dur <- which(xl_node[[1]] == "DURATION")
    row_dsw <- which(xl_node[[1]] == "DISABILITY WEIGHT")
    row_age <- which(xl_node[[1]] == "AGE")
    row_sex <- which(xl_node[[1]] == "SEX")
    
    ## setup node
    dismod_node <- list()

    ## setup node value
    dismod_node$val <- list()
    dismod_node$val$type <- xl_node[row_type, 2]  
    dismod_node$val$dist <- xl_node[row_val, 4]  
    dismod_node$val$data <- xl_node[(row_val+2):(row_dur-2), 2:6]
    colnames(dismod_node$val$data) <- xl_node[row_val+1, 2:6]
    rownames(dismod_node$val$data) <- NULL
    
    ## setup node duration
    dismod_node$dur <- list()
    dismod_node$dur$type <- "INC"
    dismod_node$dur$dist <- xl_node[row_dur, 4]  
    dismod_node$dur$data <- xl_node[(row_dur+2):(row_dsw-2), 2:6]  
    colnames(dismod_node$dur$data) <- xl_node[row_dur+1, 2:6]
    rownames(dismod_node$dur$data) <- NULL
    
    ## setup node disability weight
    dismod_node$dsw <- list()
    dismod_node$dsw$type <- "PROB"
    dismod_node$dsw$dist <- xl_node[row_dsw, 4]  
    dismod_node$dsw$data <- xl_node[(row_dsw+2):(row_age-2), 2:6]  
    colnames(dismod_node$dsw$data) <- xl_node[row_dsw+1, 2:6]
    rownames(dismod_node$dsw$data) <- NULL
    
    ## setup node age
    dismod_node$age <- list()
    dismod_node$age$type <- "AGE"
    dismod_node$age$dist <- "age"  
    dismod_node$age$data <- xl_node[(row_age+2):(row_sex-2), 2:15]  
    colnames(dismod_node$age$data) <- xl_node[row_age+1, 2:15]
    rownames(dismod_node$age$data) <- NULL
    
    ## setup node sex
    dismod_node$sex <- list()
    dismod_node$sex$type <- "SEX"
    dismod_node$sex$dist <- "sex"
    dismod_node$sex$data <- xl_node[(row_sex+2):nrow(xl_node), 2:5]  
    colnames(dismod_node$sex$data) <- xl_node[row_sex+1, 2:5]
    rownames(dismod_node$sex$data) <- NULL
    
    ## return node
    return(dismod_node)
  }

##--------------------------------------------------------------------------#
## optimizers --------------------------------------------------------------#

## Gamma distribution, for 'INC' type
optim_gamma <-
  function(par, p) {
    target <- par[c(1, 3)]
    p <- c(0, p) + (1 - p)/2
    
    f <-
      function(x, mean, p, target) {
        dev <- qgamma(p = p, shape = x, rate = x / mean)
        return(sum((dev - target) ^ 2))
      }
    
    shape <-
      optimize(f, c(0, 1000), mean = par[2], p = p, target = target)$minimum
    rate <-
      shape / par[2]
    
    return(c(shape, rate))
  }

## Beta distribution, for 'PROB' and 'SPLIT' types
optim_beta <-
  function(par, p) {
    target <- par[c(1, 3)]
    p <- c(0, p) + (1 - p)/2
    
    f <-
      function(x, mean, p, target) {
        dev <- qbeta(p = p, shape1 = x, shape2 = (x * (1 - mean))/mean)
        return(sum((dev - target) ^ 2))
      }
    
    shape1 <-
      optimize(f, c(0, 1000), mean = par[2], p = p, target = target)$minimum
    shape2 <-
      (shape1 * (1 - par[2])) / par[2]
    
    return(c(shape1, shape2))
  }

##--------------------------------------------------------------------------#
## samplers ----------------------------------------------------------------#

sim_fixed <-
  function(n, par) {
    return(rep(par, n))
  }

sim_beta <-
  function(n, par1, par2) {
    return(rbeta(n, par1, par2))
  }

sim_gamma <-
  function(n, shape, rate) {
    return(rgamma(n, shape, rate))
  }

sim_unif <-
  function(n, min, max) {
    return(runif(n, min, max))
  }

sim_pert <-
  function(n, mode, min, max) {
    bp <- prevalence::betaPERT(a = min, m = mode, b = max)
    return(rbeta(n, bp$alpha, bp$beta) * (bp$b - bp$a) + bp$a)
  }

sim_mean_ci <-
  function(n, mean, lwr, upr, type) {
    if (type == "INC") {
      pars_gamma <- optim_gamma(c(lwr, mean, upr), 0.95)
      samples <- sim_gamma(n, pars_gamma[1], pars_gamma[2])
      
    } else if (type %in% c("PROB", "SPLIT")) {
      pars_beta <- optim_beta(c(lwr, mean, upr), 0.95)
      samples <- sim_beta(n, pars_beta[1], pars_beta[2])
      
    } else {
      stop("sim_mean_ci() could not find a proper type.")
    }
    return(samples)
  }


##--------------------------------------------------------------------------#
## pre-sample input data ---------------------------------------------------#

pre_sample_input <-
  function(input, n_samples) {
    ## define seed
    ## -> deals with stratification uncertainty
    ## -> makes results reproducible
    set.seed(264)
    
    ## generate samples
    samples <-
      switch(
        input$dist,
        "Simulations" =
          mapply(
            sim_fixed,
            NA,
            n = n_samples),
        #get_simulations(n_samples, pars, type),
        "Gamma" =
          mapply(
            sim_gamma,
            as.numeric(input$data$Shape),
            as.numeric(input$data$Rate),
            n = n_samples),
        "Beta" =
          mapply(
            sim_beta,
            as.numeric(input$data$Alpha),
            as.numeric(input$data$Beta),
            n = n_samples),
        "Mean" =
          mapply(
            sim_fixed,
            as.numeric(input$data$Mean),
            n = n_samples),
        "Mean & 95%CI" =
          mapply(
            sim_mean_ci,
            as.numeric(input$data$Mean),
            as.numeric(input$data$`2.5%`),
            as.numeric(input$data$`97.5%`),
            MoreArgs = list(type = input$type),
            n = n_samples),
        "Min/Max" =
          mapply(
            sim_unif,
            as.numeric(input$data$Min),
            as.numeric(input$data$Max),
            n = n_samples),
        "Min/Mode/Max" =
          mapply(
            sim_pert,
            as.numeric(input$data$Mode),
            as.numeric(input$data$Min),
            as.numeric(input$data$Max),
            n = n_samples))
    
    ## integrate samples in dataframe
    samples_list <- do.call("c", apply(samples, 2, list))
    input_data_sim <- input$data[, c("COUNTRY", "YEAR")]
    input_data_sim$SAMPLES <- samples_list
    
    ## return samples
    return(input_data_sim)
  }

pre_sample_age <-
  function(input, n_samples) {
    ## define seed
    ## -> deals with stratification uncertainty
    ## -> makes results reproducible
    set.seed(264)
    
    ## generate samples
    age_names <- tail(names(input$data), -2)
    age_names <- age_names[age_names != "NA"]
    
    ## if sum to one, assume fixed
    if (sum(as.numeric(input$data[, age_names])) == nrow(input$data)) {
      samples <-
        apply(input$data[, age_names], 1,
              function(x) sapply(as.numeric(x), rep, each = n_samples),
              simplify = FALSE)
      samples <-
        lapply(samples, function(x) {colnames(x) <- age_names; return(x)}) 
      
    } else {
      samples <-
        apply(
          input$data[, age_names], 1, 
          function(x)
            t(rmultinom(n_samples, sum(as.numeric(x)), as.numeric(x)) /
                sum(as.numeric(x))),
          simplify = FALSE)
      samples <-
        lapply(samples, function(x) {colnames(x) <- age_names; return(x)})  
    }
    
    ## integrate samples in dataframe
    input_data_sim <- input$data[, c("COUNTRY", "YEAR")]

    for (i in seq(nrow(input_data_sim))) {
      for (j in seq_along(age_names)) {
        input_data_sim[[age_names[j]]] <-
          list(samples[[i]][, j])
      }
    }
    
    ## return samples
    return(input_data_sim)
  }

pre_sample_sex <-
  function(input, n_samples) {
    ## define seed
    ## -> deals with stratification uncertainty
    ## -> makes results reproducible
    set.seed(264)
    
    ## check for 'Both sexes'
    if (names(input$data)[3] == "Both sexes") {
      sex_names <- "Both sexes"
      samples <-
        apply(input$data[, 3, drop = FALSE], 1,
              function(x) sapply(as.numeric(x), rep, each = n_samples),
              simplify = FALSE)
      samples <-
        lapply(samples, function(x) {colnames(x) <- sex_names; return(x)}) 
      
    ## generate samples for 'Male' and 'Female'
    } else {
      sex_names <- names(input$data)[3:4]
      
      ## if sum to one, assume fixed
      if (sum(as.numeric(input$data[, 3:4])) == nrow(input$data)) {
        samples <-
          apply(input$data[, 3:4], 1,
                function(x) sapply(as.numeric(x), rep, each = n_samples),
                simplify = FALSE)
        samples <-
          lapply(samples, function(x) {colnames(x) <- sex_names; return(x)}) 
        
      } else {
        samples <-
          apply(
            input$data[, sex_names], 1, 
            function(x)
              t(rmultinom(n_samples, sum(as.numeric(x)), as.numeric(x)) /
                  sum(as.numeric(x))),
            simplify = FALSE)
        samples <-
          lapply(samples, function(x) {colnames(x) <- sex_names; return(x)})  
      }
    }
    
    ## integrate samples in dataframe
    input_data_sim <- input$data[, c("COUNTRY", "YEAR")]
    
    for (i in seq(nrow(input_data_sim))) {
      for (j in seq_along(sex_names)) {
        input_data_sim[[sex_names[j]]] <-
          list(samples[[i]][, j])
      }
    }
    
    ## return samples
    return(input_data_sim)
  }

##--------------------------------------------------------------------------#
## pre-sample node ---------------------------------------------------------#

pre_sample_node <-
  function(node, n_samples) {
    ## pre-sample only if 'dist' is specified
    node$val$samp <- 
      if (!is.na(node$val$dist)) pre_sample_input(node$val, n_samples)
    node$dur$samp <- 
      if (!is.na(node$dur$dist)) pre_sample_input(node$dur, n_samples)
    node$dsw$samp <- 
      if (!is.na(node$dsw$dist)) pre_sample_input(node$dsw, n_samples)
    node$age$samp <- 
      if (!is.na(node$age$data[1, 3])) pre_sample_age(node$age, n_samples)
    node$sex$samp <- 
      if (!is.na(node$sex$data[1, 3])) pre_sample_sex(node$sex, n_samples)
    return(node)
  }

##--------------------------------------------------------------------------#
## pre-sample dalymod list -------------------------------------------------#

pre_sample_dalymod <-
  function(dalymod, n_samples) {
    dalymod$nodes <-
      lapply(dalymod$nodes, pre_sample_node, n_samples)
    return(dalymod)
  }

##--------------------------------------------------------------------------#
## normalize samples for splits --------------------------------------------#

normalize_splits <-
  function(dalymod) {
    ## identify parents of split nodes
    parents <- unique(subset(dalymod$dismod, TYPE == "SPLIT")$PARENT)
    
    ## normalize for each split group
    for (i in seq_along(parents)) {
      ## get split nodes
      nodes <- subset(dalymod$dismod, PARENT == parents[i])$NODE
      
      ## compile samples in nested list
      samp_list <- lapply(nodes, function(x) dalymod$nodes[[x]]$val$samp$SAMPLES)
      n_inner <- length(samp_list[[1]]) # number of country-years
      n_outer <- length(samp_list)      # number of nodes
      
      ## loop over samples and normalize
      ## replace original with normalized samples
      for (j in seq(n_inner)) {
        samp_list_sum <- rowSums(sapply(samp_list, function(x) x[[j]]))
        
        for (k in seq(n_outer)) {
          dalymod$nodes[[nodes[k]]]$val$samp$SAMPLES[[j]] <-
            samp_list[[k]][[j]] / samp_list_sum
        }
      }
    }
    
    ## return normalized dalymod
    return(dalymod)
  }

##--------------------------------------------------------------------------#
## main wrapper functions --------------------------------------------------#

dalymod <-
  function(file, n_samples) {
    ## import dalymod excel
    dalymod <- import_dalymod(file)
    
    ## pre-sample nodes
    dalymod <- pre_sample_dalymod(dalymod, n_samples)
    
    ## normalize splits
    dalymod <- normalize_splits(dalymod)
    
    ## multiply nodes to get incidence per terminal node
    dalymod <- multiply_dalymod(dalymod)
    
    ## return updated dalymod
    return(dalymod)
  }

file <-
  "C:/Users/BrDe394/Dropbox/#FERG2/#CTF/daly-calculation/ferg2-daly-anthrax.xlsx"
dalymods <- dalymod(file, 5)

## set year
year <- 2010

## aggregate population by age groups
pop <- subset(FERG2:::pop, YEAR == year)
pop$AGE2 <- cut(pop$AGE, c(seq(0, 85, 5), Inf), right = FALSE)
pop2 <- aggregate(POP ~ AGE2 + SEX + ISO3, pop, sum)
names(pop2)[names(pop2) == "AGE2"] <- "AGE"
str(pop2)

dalycalc <-
  function(dalymod, year, pop) {
    
  }


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
  function(agesex_agg, country) {
    agesex_agg_pop <-
      expand.grid(
        AGE = split_age_string(agesex_agg$AGE),
        SEX = split_sex_string(agesex_agg$SEX))
    agesex_agg_pop <-
      merge(agesex_agg_pop, subset(pop2, ISO3 == country))
    agesex_agg_pop$W <-
      agesex_agg_pop$POP / sum(agesex_agg_pop$POP)
    agesex_agg_pop$INC_NR <-
      lapply(agesex_agg_pop$W, function(x) x * agesex_agg$INC_NR)
    
    agesex_agg_pop$W <- NULL
    agesex_agg_pop$ISO3 <- NULL
    
    apply(agesex_agg_pop, 1, as.list)
  }

split_agesex_all <-
  function(agesex_agg_all, country) {
    agesex_all <- sapply(agesex_agg_all, split_agesex, "ZWE")
    agesex_all <- unlist(agesex_all, recursive = FALSE)
    return(agesex_all)
  }



node <- dalymods$nodes$split_cutaneous_mild

## TO DO
## .. add contribution/inc to nodes
## .. calculate YLD, YLL ifo contrib
## .. .. calculate rle per age/sex
## .. apply dalycalc_node to all nodes > dalycalc
## .. calculate summaries across nodes > sum yld, yll, inc by country
## .. .. summary ~ age, sex, ..
## .. calculate aggregates across countries

dalycalc_node <-
  function(node, year, pop) {
    ## aggregate pop by country
    pop_agg <- aggregate(POP ~ ISO3, pop, sum)

    ## expand countries if needed
    if (node$inc$COUNTRY == "ALL") {
      node_inc <- cbind(pop_agg$ISO3, node$inc)
      node_inc$COUNTRY <- NULL
      names(node_inc)[names(node_inc) == "pop_agg$ISO3"] <- "COUNTRY"
    }

    ## add year if needed
    if (node$inc$YEAR == "ALL") {
      node_inc$YEAR <- year
    }
    
    ## merge incidence and population
    node_inc <- merge(node_inc, pop_agg, by.x = "COUNTRY", by.y = "ISO3")
    
    ## multiply inc with pop
    node_inc$INC_NR <-
      apply(node_inc, 1,
            function(x) x[["POP"]] * x[["SAMPLES"]], simplify = FALSE)
    
    ## setup age split
    age_split <- node$age$samp
    age_split$YEAR <- NULL
    if (age_split$COUNTRY == "ALL") {
      age_split <- cbind(pop_agg$ISO3, age_split)
      age_split$COUNTRY <- NULL
      names(age_split)[names(age_split) == "pop_agg$ISO3"] <- "COUNTRY"
    }
    rownames(age_split) <- age_split$COUNTRY
    age_split$COUNTRY <- NULL
    age_names <- names(age_split)
    
    ## setup sex split
    sex_split <- node$sex$samp
    sex_split$YEAR <- NULL
    if (sex_split$COUNTRY == "ALL") {
      sex_split <- cbind(pop_agg$ISO3, sex_split)
      sex_split$COUNTRY <- NULL
      names(sex_split)[names(sex_split) == "pop_agg$ISO3"] <- "COUNTRY"
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
      function(x) split_agesex_all(dalycalc_agg[[x]], x))
    names(dalycalc_all) <- names(dalycalc_agg)

    dalycalc_all[[1]]
    
    ## calculate YLD
    if (node$contribution == "YLD") {
      for (i in seq_along(dalycalc_all)) {
        for (j in seq_along(dalycalc_all[[i]])){
          dalycalc_all[[i]][[j]]$YLD_NR <-
            dalycalc_all[[i]][[j]]$INC_NR * 
            node$dur$samp$SAMPLES *
            node$dsw$samp$SAMPLES
        }
      }
    }

    ## calculate YLL
    if (node$contribution == "YLL") {
      for (i in seq_along(dalycalc_all)) {
        for (j in seq_along(dalycalc_all[[i]])){
          dalycalc_all[[i]][[j]]$YLL_NR <-
            dalycalc_all[[i]][[j]]$INC_NR * 
            node$dur$samp$SAMPLES
        }
      }
    }
    
    ## return dalycalc node
    return(dalycalc_all)
}

             
    ## ISO > 18*2 > AGE,SEX,POP,INC_NR,INC_RT,YLD_NR,YLD_RT,YLL_NR,YLL_RT

    ## calculate YLD
    ## calculate YLL
    ## calculate INC
    

multiply_dalymod <-
  function(dalymod) {
    ## extract all nodes
    nodes <- dalymod$dismod$NODE
    
    ## per node, multiply across tree
    for (i in seq_along(nodes)) {
      dalymod$nodes[[nodes[i]]]$inc <- multiply_nodes(dalymod, nodes[i])
    }
    
    ## return dalymod
    return(dalymod)
  }

multiply_nodes <-
  function(dalymod, node) {
    ## get tree
    tree <- get_tree(dalymod$dismod, node)
    
    ## compile samples in nested list
    samp_list <- lapply(tree, function(x) dalymod$nodes[[x]]$val$samp$SAMPLES)
    n_inner <- length(samp_list[[1]]) # number of country-years

    ## prepare output
    inc <- dalymod$nodes[[node]]$val$samp[, 1:2]

    ## loop over samples and calculate product
    for (j in seq(n_inner)) {
      samp_list_prod <- apply(sapply(samp_list, function(x) x[[j]]), 1, prod)
      inc$SAMPLES[j] <- list(samp_list_prod)
    }
    
    ## return inc
    return(inc)
  }

get_tree <-
  function(dismod, node) {
    tree <- parent <- node
    while(parent != "NA") {
      parent <- subset(dismod, NODE == parent)$PARENT
      tree <- c(tree, parent)
    }
    return(head(tree, -1))
  }




dalycalc

n_samples <- 1e5

expand.grid(
  COUNTRY = c("BEL", "ETH"),
  YEAR = 2000:2020,
  AGE = 1:18,
  SEX = 1:2)

INC/DUR/DSW/YLD
INC/RLE/YLL

dalymod

dalycalc <- list()
dalycalc$inc <- vector("list", length(nodes))
dalycalc$yld <- vector("list", length(nodes))
dalycalc$yll <- vector("list", length(nodes))

names(dalycalc$inc) <- nodes
dalycalc$yld$split_cutaneous_mild <-
  inc

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

