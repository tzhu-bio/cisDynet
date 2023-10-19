## these functions were adapted from the cicero package

#' Generate window for a given bed file.
#'
#' @param win
#' @param genomic_coords
#'
#' @return
#' @export
#'
#' @examples
generate_windows <- function(win, genomic_coords) {
  if(!is(genomic_coords, "data.frame")) {
    chr_maxes <- read.table(genomic_coords)
  } else {
    chr_maxes <- genomic_coords
  }
  names(chr_maxes) <- c("V1", "V2")
  win_ranges <- plyr::ddply(chr_maxes, plyr::.(V1), function(x) {
    r <- seq(from = 1, to = x$V2[1], by = win/2)
    l <- r + win - 1
    data.frame(start = r, end = l)
  })
  gr <- GenomicRanges::GRanges(win_ranges$V1,
                               ranges=IRanges::IRanges(win_ranges$start,
                                                       win_ranges$end))
  return(gr)
}

get_vals <- function(dat){
  dat[, !(colnames(dat) %in% c('chr','bp1','bp2','summit','mean_bp'))]
}

get_genomic_range <- function(grs, indata, win) {
  end1 <- as.numeric(as.character(GenomicRanges::end(grs[win])))
  end2 <- as.numeric(as.character(GenomicRanges::start(grs[win])))
  win_range <- indata[(indata$bp1 < end1 &
                         indata$bp1 > end2) |
                        (indata$bp2 < end1 &
                           indata$bp2 > end2), ]
  # win_range <-
  #   win_range[as.character(win_range$chr) ==
  #               gsub("chr", "",
  #                    as.character(GenomicRanges::seqnames(grs[win]))),]
  win_range <-
    win_range[as.character(win_range$chr) == as.character(GenomicRanges::seqnames(grs[win])),]
  win_range$mean_bp <-
    (as.numeric(as.character(win_range$bp1)) +
       as.numeric(as.character(win_range$bp2)))/2

  return(win_range)
}



get_rho_mat <- function(dist_matrix, distance_parameter, s) {
  xmin <- 1000
  out <- (1-(xmin/dist_matrix)^s) * distance_parameter
  out[!is.finite(out)] <- 0
  out[out < 0] <- 0
  return(out)
}

find_distance_parameter <- function(dist_mat,
                                    gene_range,
                                    maxit,
                                    null_rho,
                                    s,
                                    distance_constraint,
                                    distance_parameter_convergence) {
  if (sum(dist_mat > distance_constraint)/2 < 1) {
    #warning("No long edges")
    return("No long edges")
  }

  found <- FALSE
  starting_max <- 2
  distance_parameter <- 2
  distance_parameter_max <- 2
  distance_parameter_min <- 0
  it <- 0
  while(found != TRUE & it < maxit) {
    vals <- get_vals(gene_range)
    cov_mat <- cov(t(vals))
    diag(cov_mat) <- diag(cov_mat) + 1e-4

    rho <- get_rho_mat(dist_mat, distance_parameter, s)

    GL <- glasso::glasso(cov_mat, rho)
    big_entries <- sum(dist_mat > distance_constraint)

    if (((sum(GL$wi[dist_mat > distance_constraint] != 0)/big_entries) > 0.05) |
        (sum(GL$wi == 0)/(nrow(GL$wi)^2) < 0.2 ) ) {
      longs_zero <- FALSE
    } else {
      longs_zero <- TRUE
    }

    if (longs_zero != TRUE | (distance_parameter == 0)) {
      distance_parameter_min <- distance_parameter
    } else {
      distance_parameter_max <- distance_parameter
    }
    new_distance_parameter <- (distance_parameter_min +
                                 distance_parameter_max)/2

    if(new_distance_parameter == starting_max) {
      new_distance_parameter <- 2 * starting_max
      starting_max <- new_distance_parameter
    }

    if (distance_parameter_convergence > abs(distance_parameter -
                                             new_distance_parameter)) {
      found <- TRUE
    } else {
      distance_parameter <- new_distance_parameter
    }
    it <- it + 1
  }
  if (maxit == it) warning("maximum iterations hit")
  return(distance_parameter)
}


#' Calculate distance for matrix.
#'
#' @param gene_range
#'
#' @return
#' @export
#'
#' @examples
calc_dist_matrix <- function(gene_range) {
  dist_mat <- as.matrix(dist(gene_range$mean_bp))
  row.names(dist_mat) <- colnames(dist_mat) <- row.names(gene_range)
  return(dist_mat)
}


#' Title
#'
#' @param indata
#' @param window
#' @param maxit
#' @param s
#' @param sample_num
#' @param distance_constraint
#' @param distance_parameter_convergence
#' @param max_elements
#' @param genomic_coords
#'
#' @return
#' @export
#'
#' @examples
estimate_distance_parameter <- function(indata,
                                        window=window,
                                        maxit=100,
                                        s=0.75,
                                        sample_num = 100,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        max_elements = 200,
                                        genomic_coords = NULL) {

  assertthat::assert_that(is(genomic_coords, 'data.frame'))

  grs <- generate_windows(window, genomic_coords)

  distance_parameters <- list()
  distance_parameters_calced <- 0
  it <- 0

  while(sample_num > distance_parameters_calced & it < 3 * sample_num) {
    it <- it + 1
    win <- sample(seq_len(length(grs)), 1)
    GL <- "Error"
    win_range <- get_genomic_range(grs, indata, win)

    if (nrow(win_range)<=1) {
      next()
    }
    if (nrow(win_range) > max_elements) {
      next()
    }

    dist_matrix <- calc_dist_matrix(win_range)

    distance_parameter <- find_distance_parameter(dist_matrix,
                                                  win_range,
                                                  maxit = maxit,
                                                  null_rho = 0,
                                                  s,
                                                  distance_constraint = distance_constraint,
                                                  distance_parameter_convergence = distance_parameter_convergence)

    if (!is(distance_parameter, "numeric")) next()
    distance_parameters = c(distance_parameters, distance_parameter)
    distance_parameters_calced <- distance_parameters_calced + 1
  }

  #if(length(distance_parameters) < sample_num)
    #warning(paste("Could not calculate sample_num distance_parameters - see",
    #              "documentation details", collapse = " "))
  if(length(distance_parameters) == 0)
    stop("No distance_parameters calculated")

  unlist(distance_parameters)
}

#' Title
#'
#' @param indata
#' @param distance_parameter
#' @param s
#' @param window
#' @param max_elements
#' @param genomic_coords
#'
#' @return
#' @export
#'
#' @examples
generate_cicero_models <- function(indata,
                                   distance_parameter,
                                   s = 0.75,
                                   window = 500000,
                                   max_elements = 200,
                                   genomic_coords = NULL) {

  assertthat::assert_that(is(genomic_coords, 'data.frame'))

  grs <- generate_windows(window, genomic_coords)

  outlist <- parallel::mclapply(seq_len(length(grs)), mc.cores = 1, function(win) {
    GL <- "Error"

    win_range <- get_genomic_range(grs, indata, win)

    if (nrow(win_range)<=1) {
      #warning("Zero or one element in range")
      return("Zero or one element in range")
    }
    if (nrow(win_range) > max_elements) {
      #warning("Too many elements in range")
      return("Too many elements in range")
    }

    dist_matrix <- calc_dist_matrix(win_range)

    rho_mat <- get_rho_mat(dist_matrix, distance_parameter, s)

    vals <- get_vals(win_range)
    cov_mat <- cov(t(vals))
    diag(cov_mat) <- diag(cov_mat) + 1e-4

    GL <- glasso::glasso(cov_mat, rho_mat)
    colnames(GL$w) <- row.names(GL$w) <- row.names(vals)
    colnames(GL$wi) <- row.names(GL$wi) <- row.names(vals)
    return(GL)
  })
  names_df <- as.data.frame(grs)
  names(outlist) <- paste(names_df$seqnames,
                          names_df$start,
                          names_df$end, sep="_")

  #FIXME add warning about how many regions removed due to too many elements
  outlist
}

#' Title
#'
#' @param values
#'
#' @return
#' @export
#'
#' @examples
reconcile <- function(values) {
  if (length(values) == 1) return(values)
  if (sum(values >= 0) == length(values)) return(mean(values))
  if (sum(values <= 0) == length(values)) return(mean(values))
  if (sum(values == 0) == length(values)) return(0)
  return(NA_real_)
}

#' Title
#'
#' @param cicero_model_list
#' @param silent
#'
#' @return
#' @export
#'
#' @examples
assemble_connections <- function(cicero_model_list, silent = FALSE) {
  types <- vapply(cicero_model_list, FUN=class, FUN.VALUE="character")
  char_hbn <- cicero_model_list[types=="character"]
  gl_only <- cicero_model_list[types=="list"]
  if(!silent) {
    print(paste("Successful cicero models: ", length(gl_only)))
    print("Other models: ")
    print(table(unlist(char_hbn)))
    print(paste("Models with errors: ", sum(is.null(cicero_model_list))))
  }

  cors <- lapply(gl_only, function(gl)  {
    cors <- stats::cov2cor(gl$w)
    data.table::melt(data.table::as.data.table(cors, keep.rownames=TRUE),
                     measure=patterns("[0-9]"))
  })

  cors <- data.table::rbindlist(cors)
  names(cors) <- c("Var1", "Var2", "value")
  data.table::setkey(cors, "Var1", "Var2")

  cors_rec <- cors #as.data.frame(cors[,list(mean_coaccess = reconcile(value)), by="Var1,Var2"])

  names(cors_rec) <- c("Peak1", "Peak2", "coaccess")
  cors_rec <- cors_rec[cors_rec$Peak1 != cors_rec$Peak2,]
  return(cors_rec)
}

#' Title
#'
#' @param indata
#' @param genomic_coords
#' @param window
#' @param silent
#' @param sample_num
#'
#' @return
#' @export
#'
#' @examples
run_cicero <- function(indata,
                       genomic_coords,
                       window = 500000,
                       silent=FALSE,
                       sample_num = 100) {
  # Check input
  assertthat::assert_that(is(genomic_coords, 'data.frame'))

  #if (!silent) logfile("Starting ...")
  #if (!silent) logfile("Calculating distance_parameter value")
  distance_parameters <- estimate_distance_parameter(indata, window=window,
                                                     maxit=100, sample_num = sample_num,
                                                     distance_constraint = 250000,
                                                     distance_parameter_convergence = 1e-22,
                                                     genomic_coords = genomic_coords)

  mean_distance_parameter <- mean(unlist(distance_parameters))

  #if (!silent) logfile("Running models")
  cicero_out <-
    generate_cicero_models(indata,
                           distance_parameter = mean_distance_parameter,
                           window = window,
                           genomic_coords = genomic_coords)
  #return(cicero_out)
  #logfile("Assembling connections")
  all_cons <- assemble_connections(cicero_out, silent=TRUE)

  #logfile("Done")
  return(all_cons)
}
