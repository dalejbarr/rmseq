#' Basic functions for trial sequencing
#'
#' These functions reflect alternative ways of generating constrained
#' random trial sequences.
#'
#' @name sequencing-functions
#' @examples
#' # TODO
NULL

pad_id <- function(n, prefix, factor = TRUE) {
  fmt_str <- paste0(prefix, "%0", ceiling(log(n + 1L, 10L)), "d")
  result <- sprintf(fmt_str, seq_len(n))
  if (factor) {
    factor(result)
  } else {
    result
  }
}  

check_names <- function(x) {
  noun <- if (is.list(x)) "list" else "vector"
  if (is.null(names(x))) {
    stop("'labels' must be a named ", noun)
  } else if (any(names(x) == "") && !is.list(x)) {
    stop("All elements of 'labels' must be named.")
  }

  if (length(unique(names(x))) != length(names(x))) {
    stop("Elements of 'labels' must have unique names.")
  }
}

generate_factor_labels <- function(x, nx) {
  if (is.numeric(x)) {
    if (x < 2L) {
      stop("Number of factor levels must be greater than 1.")
    }
    res <- paste0(nx, seq_len(x))
  } else {
    res <- x
  }
  xf <- factor(res)
  xf
}

## process user-supplied input to 'labels' argument in sequencing function
## and check them
is_factorial_input <- function(labels) {
  is_factorial <- FALSE

  if (is.list(labels)) {
    check_names(labels)
    is_factorial <- TRUE
  } else {
    ## it's a vector
    if (all(sapply(labels, is.numeric)) &&
        all(sapply(labels, length) == 1L)) {
      check_names(labels)
      is_factorial <- TRUE
    }
  }

  is_factorial
}

normalize_factorial_input <- function(labels) {
  mapply(generate_factor_labels, labels, names(labels),
         SIMPLIFY=FALSE)
}

## process user-supplied input to 'labels' argument in sequencing function
## and return a set of levels
get_levels <- function(labels) {
  if (is_factorial_input(labels)) {
    faclist <- normalize_factorial_input(labels)
    suppressWarnings(levels(interaction(faclist, sep = ":")))
  } else {
    labels
  }
}

#' @rdname sequencing-functions
#'
#' @param labels A vector containing levels of a nominal variable or a
#'   list specifying multiple variables in a factorial design (see
#'   'Details').
#'
#' @param n_buckets Number of 'buckets' in which repeated labels are
#'   to be organized.
#'
#' @param n_reps Number of repetitions of each label per bucket.
#'
#' @param edge_streaks Whether to allow streaks at bucket edges.
#'
#' @param repeats Whether each bucket repeats a single random sequence.
#'
#' @details TODO.
#'
#' @return A vector of length `n_buckets * length(labels) * n_reps`
#'   with a constrained sequence of factor labels (or cell labels for
#'   factorial designs).
#'
#' @export
seq_n <- function(labels, n_buckets, n_reps = 1,
                  edge_streaks = FALSE, repeats = FALSE) {
  
  .seq_n(get_levels(labels), n_buckets, n_reps, edge_streaks, repeats)
}

#' @rdname sequencing-functions
#' @inheritParams seq_n
#' @export
com_n <- function(labels, n_buckets, n_reps = 1,
                  edge_streaks = FALSE, repeats = FALSE) {

  .com_n(get_levels(labels), n_buckets, n_reps, edge_streaks, repeats)
}

#' @rdname sequencing-functions
#' @inheritParams seq_n
#' @param n_subj Number of subjects.
#' @param comp_sets Whether to generate complement sets.
#' @export
plan_n <- function(n_subj,
                   labels, n_buckets, n_reps = 1,
                   edge_streaks = FALSE, repeats = FALSE,
                   comp_sets = FALSE) {
  
  design <- if (is_factorial_input(labels)) {
              normalize_factorial_input(labels)
            } else {
              list(cond = generate_factor_labels(labels, "cond"))
            }

  combos <- expand.grid(lapply(rev(design), levels))
  design_table <- combos[rev(seq_along(combos))]

  ## TODO: check id, trial, cseq
  if (any(names(design_table) == ".key")) {
    stop("Variables cannot be named '.key'")
  }
  design_table[[".key"]] <- seq_len(nrow(design_table))

  study_plan <-
    if (comp_sets) {
      if ((n_subj %% nrow(design_table)) > 0L) {
        stop("When 'comp_sets' is TRUE, 'n_subj' must be a multiple of\n",
             "  the number of labels (or cells in a factorial design), which\n",
             "  is ", nrow(design_table), " for the current data.")
      }
      plan_lists <- replicate(n_subj / nrow(design_table),
                              .com_n_ix(seq_len(nrow(design_table)),
                                        n_buckets, n_reps,
                                        edge_streaks, repeats),
                              simplify = FALSE)
      do.call("c", plan_lists)
    } else {
      replicate(n_subj, .seq_n_ix(seq_len(nrow(design_table)),
                                  n_buckets, n_reps,
                                  edge_streaks, repeats),
                simplify = FALSE)
    }


  ntrials <- length(study_plan[[1]]) # number trials per subject
  ncompsets <- n_subj / nrow(design_table)
  stbl <- data.frame(id = rep(pad_id(n_subj, "S"), each = ntrials))
  if (comp_sets) {
    stbl[["cset"]] <- paste(rep(pad_id(ncompsets, "C"),
                                each = ntrials * nrow(design_table)),
                            rep(rep(seq_len(nrow(design_table)), each = ntrials),
                                times = ncompsets),
                            sep = ".")
  }
  stbl[["trial"]] <- rep(seq_along(study_plan[[1]]), times = n_subj)
  stbl[[".key"]] <- unlist(study_plan)

  result <- merge(stbl, design_table, by = ".key", sort = FALSE)
  result_ord <- result[order(result$id, result$trial), ]
  rownames(result_ord) <- NULL
  result_ord[setdiff(names(result_ord), ".key")]
}
