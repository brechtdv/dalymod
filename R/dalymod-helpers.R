###=========================================================================#
### dalymod helpers 
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
###---------| get_samples ... actual sampler, get simulations from file
###-----------| expit ....... inverse logit
###---------| sim_fixed ..... actual sampler, fixed
###---------| sim_beta ...... actual sampler, beta
###---------| sim_gamma ..... actual sampler, gamma
###---------| sim_unif ...... actual sampler, uniform
###---------| sim_pert ...... actual sampler, beta-pert
###---------| sim_mean_ci ... actual sampler, mean & 95%CI
###---------| optim_gamma ... find gamma pars through optimization
###---------| optim_beta .... find beta pars through optimization
###---------| betaPert ...... calculate beta pars for PERT
###-- normalize_splits ...... normalize samples for splits
###-- multiply_dalymod ...... multiply all nodes in dalymod object
###---| multiply_nodes ...... multiply node in dalymod object
###-----| get_tree .......... get parent nodes from disease model
###-- dalymod ............... main wrapper function

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
    #names(xl_nodes) <- tolower(nodes)
    names(xl_nodes) <- nodes
    for (i in seq_along(nodes)) {
      xl_nodes[[i]] <- import_node(file, nodes[i])
    }
    
    ## return import
    return(list(dismod = xl_dismod, nodes = xl_nodes))
  }

import_node <-
  function(file, node) {
    ## read Excel node
    xl_node <-
      read_excel_base(file, node, col_names = FALSE, col_types = "text")
    
    ## define row numbers
    row_type <- which(xl_node[[1]] == "TYPE")
    row_den <- which(xl_node[[1]] == "DENOMINATOR")
    row_cnt <- which(xl_node[[1]] == "CONTRIBUTION")
    row_inc <- which(xl_node[[1]] == "INCIDENCE")
    row_val <- which(xl_node[[1]] == "VALUE")
    row_dur <- which(xl_node[[1]] == "DURATION")
    row_dsw <- which(xl_node[[1]] == "DISABILITY WEIGHT")
    row_age <- which(xl_node[[1]] == "AGE")
    row_sex <- which(xl_node[[1]] == "SEX")
    
    ## setup node
    dismod_node <- list()

    ## setup node settings
    dismod_node$set$contribution <- xl_node[row_cnt, 2]
    dismod_node$set$incidence <- xl_node[row_inc, 2]

    ## setup node value
    dismod_node$val <- list()
    dismod_node$val$type <- xl_node[row_type, 2]  
    dismod_node$val$dist <- xl_node[row_val, 4]  
    dismod_node$val$data <- xl_node[(row_val+2):(row_dur-2), 2:6]
    dismod_node$val$dnmn <- as.numeric(xl_node[row_den, 2])
    colnames(dismod_node$val$data) <- xl_node[row_val+1, 2:6]
    rownames(dismod_node$val$data) <- NULL
    
    ## setup node duration
    dismod_node$dur <- list()
    dismod_node$dur$type <- "INC"
    dismod_node$dur$dist <- xl_node[row_dur, 4]  
    dismod_node$dur$data <- xl_node[(row_dur+2):(row_dsw-2), 2:6]
    dismod_node$dur$dnmn <- 1
    colnames(dismod_node$dur$data) <- xl_node[row_dur+1, 2:6]
    rownames(dismod_node$dur$data) <- NULL
    
    ## setup node disability weight
    dismod_node$dsw <- list()
    dismod_node$dsw$type <- "PROB"
    dismod_node$dsw$dist <- xl_node[row_dsw, 4]  
    dismod_node$dsw$data <- xl_node[(row_dsw+2):(row_age-2), 2:6]
    dismod_node$dsw$dnmn <- 1
    colnames(dismod_node$dsw$data) <- xl_node[row_dsw+1, 2:6]
    rownames(dismod_node$dsw$data) <- NULL
    
    ## setup node age
    dismod_node$age <- list()
    dismod_node$age$type <- "AGE"
    dismod_node$age$dist <- "age"  
    dismod_node$age$data <- xl_node[(row_age+2):(row_sex-2), 2:21]  
    colnames(dismod_node$age$data) <- xl_node[row_age+1, 2:21]
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

## calculate beta pars for PERT
betaPERT <-
function (a, m, b, k = 4, method = c("classic", "vose")) {
  if (!exists("a")) 
      stop("'a' is missing")
  if (!exists("m")) 
      stop("'m' is missing")
  if (!exists("b")) 
      stop("'b' is missing")
  if (!is.numeric(a)) 
      stop("'a' must be a numeric value")
  if (!is.numeric(m)) 
      stop("'m' must be a numeric value")
  if (!is.numeric(b)) 
      stop("'b' must be a numeric value")
  if (!exists("method")) 
      stop("'method' is missing")
  method <- match.arg(method)
  if (method == "classic") {
      if (!exists("k")) 
          stop("'k' is missing")
      if (!is.numeric(k)) 
          stop("'k' must be a numeric value")
      mu <- (a + k * m + b)/(k + 2)
      sdev <- (b - a)/(k + 2)
      alpha <- ((mu - a)/(b - a)) * (((mu - a) * (b - mu)/(sdev^2)) - 
          1)
      beta <- alpha * (b - mu)/(mu - a)
  }
  if (method == "vose") {
      if (!exists("k")) 
          stop("'k' is missing")
      if (!is.numeric(k)) 
          stop("'k' must be a numeric value")
      mu <- (a + k * m + b)/(k + 2)
      alpha <- ifelse(mu == m, 1 + k/2, ((mu - a) * (2 * m - 
          a - b))/((m - mu) * (b - a)))
      beta <- alpha * (b - mu)/(mu - a)
  }
  out <- list(alpha = alpha, beta = beta, a = a, m = m, b = b, 
      method = method)
  class(out) <- "betaPERT"
  return(out)
}


##--------------------------------------------------------------------------#
## samplers ----------------------------------------------------------------#

## make several functions?
## get_samples_sci
## get_samples_who
## get_samples_gbd

expit <-
  function(x) {
    exp(x) / (1 + exp(x))
  }

get_samples <-
function(n_samples, file, transformation, denominator) {
  ## import file
  f <- paste0("../ESTIMATES/", file)
  if (!file.exists(f))
    stop("File ", sQuote(file), " not found in 'ESTIMATES' folder.")
  sim <- readRDS(f)

  ## add year 'ALL' if year is missing
  if (is.null(sim$YEAR)) sim$YEAR <- "ALL"
  
  ## compile output containing country/year/samples
  sim$SAMPLES <-
    sim$SIM[, sample(seq(ncol(sim$SIM)), n_samples, replace = TRUE), drop = F]
  sim$SAMPLES[sim$SAMPLES == 0] <- NA  # manually set to NA to maintain zero's
  sim$SAMPLES <-
    switch(
      transformation,
      log = exp(sim$SAMPLES),
      logit = expit(sim$SAMPLES))
  sim$SAMPLES <- sim$SAMPLES / as.numeric(denominator)
  sim$SAMPLES[is.na(sim$SAMPLES)] <- 0  # manually set back to 0
  sim$SAMPLES <- apply(sim$SAMPLES, 1, c, simplify = FALSE)
  
  ## expand output as needed
  sim_out <-
    expand.grid(
      ISO3 = FERG2:::countries$ISO3,
      YEAR = unique(sim$YEAR))
  sim_out <-
    base::merge(sim_out, FERG2:::countries[, c("ISO3", "REG2", "SUB2")])
  names(sim_out)[names(sim_out) == "ISO3"] <- "COUNTRY"
  
  ## integrate samples
  if (!is.null(sim$COUNTRY)) {
    sim_out <-
      base::merge(sim_out, sim[, c("COUNTRY", "YEAR", "SAMPLES")])
    
  } else if (!is.null(sim$SUB2)) {
    sim_out <-
      base::merge(sim_out, sim[, c("SUB2", "YEAR", "SAMPLES")])
    
  } else if (!is.null(sim$REG2)) {
    sim_out <-
      base::merge(sim_out, sim[, c("REG2", "YEAR", "SAMPLES")])
    
  } else {
    sim_out <-
      base::merge(sim_out, sim[, c("YEAR", "SAMPLES")])
  }

  sim_out <- sim_out[, c("COUNTRY", "YEAR", "SAMPLES")]

  ## return output
  return(sim_out)
}

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
    bp <- betaPERT(a = min, m = mode, b = max)
    return(rbeta(n, bp$alpha, bp$beta) * (bp$b - bp$a) + bp$a)
  }

sim_mean_ci <-
  function(n, mean, lwr, upr, type) {
    if (type %in% c("INC", "RATIO")) {
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

    ## take samples from file
    if (input$dist == "Simulations")
      return(get_samples(
        n_samples,
        input$data$`File name`,
        input$data$Transformation,
        input$dnmn))
    
    ## lifelong duration
    if (input$dist == "Lifelong")
      return(
        data.frame(
          COUNTRY = "ALL",
          YEAR = "ALL",
          LIFELONG = TRUE))

    ## generate samples
    samples <-
      switch(
        input$dist,
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
        "Mode/Min/Max" =
          mapply(
            sim_pert,
            as.numeric(input$data$Mode),
            as.numeric(input$data$Min),
            as.numeric(input$data$Max),
            n = n_samples))
    
    ## replace any 'NaN' by 0
    if (any(is.nan(samples))) {
      samples[is.nan(samples)] <- 0
      warning("NaN values replaced by 0")
    }
    
    ## divide by denominator
    samples <- samples / input$dnmn
    
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
    age_names <- age_names[!is.na(age_names)]    
    age_names <- age_names[age_names != "COMMENT"]
    
    ## if sum to one, assume fixed
    if (isTRUE(all.equal(
          sum(as.numeric(as.matrix(input$data[, age_names]))),
          nrow(input$data)))) {
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
      if (isTRUE(all.equal(
            sum(as.numeric(as.matrix(input$data[, 3:4]))),
            nrow(input$data)))) {
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
  function(.dalymod, n_samples) {
    .dalymod$nodes <-
      lapply(.dalymod$nodes, pre_sample_node, n_samples)
    return(.dalymod)
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
## multiply values across tree ---------------------------------------------#

multiply_dalymod <-
  function(.dalymod) {
    ## extract all nodes
    nodes <- .dalymod$dismod$NODE
    
    ## per node, multiply across tree
    for (i in seq_along(nodes)) {
      .dalymod$nodes[[nodes[i]]]$inc <- multiply_nodes(.dalymod, nodes[i])
    }
    
    ## return dalymod
    return(.dalymod)
  }

multiply_nodes <-
  function(.dalymod, node) {
    ## get tree
    tree <- rev(get_tree(.dalymod$dismod, node))
      
    ## prepare samples dataframe
    samp_df <- .dalymod$nodes[[1]]$val$samp
    if (all(samp_df$COUNTRY == "ALL")) samp_df$COUNTRY <- NULL
    if (all(samp_df$YEAR == "ALL")) samp_df$YEAR <- NULL
    
    ## merge in possible other samples
    for (i in seq_along(tree)[-1]) {
      samp_df2 <- .dalymod$nodes[[tree[i]]]$val$samp
      if (all(samp_df2$COUNTRY == "ALL")) samp_df2$COUNTRY <- NULL
      if (all(samp_df2$YEAR == "ALL")) samp_df2$YEAR <- NULL
      names(samp_df2)[names(samp_df2) == "SAMPLES"] <- paste0("SAMPLES", i)
      samp_df <- base::merge(samp_df, samp_df2)
    }
    
    ## add COUNTRY and YEAR if needed
    if (is.null(samp_df$COUNTRY))
      samp_df$COUNTRY <- "ALL"
    if (is.null(samp_df$YEAR))
      samp_df$YEAR <- "ALL"
    
    ## prepare output
    inc <- samp_df[, c("COUNTRY", "YEAR")]
    n_inner <- nrow(inc)
    
    ## loop over samples and calculate product
    for (j in seq(n_inner)) {
      samp_list_prod <-
        apply(
          sapply(
            samp_df[, grepl("SAMPLES", names(samp_df)), drop = FALSE],
            function(x) x[[j]]),
          1, prod)
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

##--------------------------------------------------------------------------#
## main wrapper functions --------------------------------------------------#

dalymod <-
  function(file, n_samples, verbose = TRUE) {
    ## set cli 'cli.default_handler' option
    if (verbose) {
      options(cli.default_handler = NULL)
    } else {
      options(cli.default_handler = function(...) { })
    }
    
    ## import dalymod excel
    cli_progress_step("Importing file {.file {file}}", spinner = TRUE)
    .dalymod <- import_dalymod(file)

    ## pre-sample nodes
    cli_progress_step("Pre-sampling nodes", spinner = TRUE)
    .dalymod <- pre_sample_dalymod(.dalymod, n_samples)

    ## normalize splits
    cli_progress_step("Normalizing samples", spinner = TRUE)
    .dalymod <- normalize_splits(.dalymod)

    ## multiply nodes to get incidence per terminal node
    cli_progress_step("Muliplying samples across nodes", spinner = TRUE)
    .dalymod <- multiply_dalymod(.dalymod)
    
    ## return updated dalymod
    return(.dalymod)
  }
