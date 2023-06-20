#!/usr/bin/env Rscript

### import packages
invisible(suppressMessages(lapply(
  c("ggplot2", "ggseqlogo", "scales", "cowplot", "optparse"),
  require, character.only = TRUE)))

### themes
no_axes    <- list(theme(line            = element_blank(),
                         rect            = element_blank(),
                         text            = element_blank(),
                         axis.line       = element_blank(),
                         axis.text.x     = element_blank(),
                         axis.ticks      = element_blank(),
                         legend.position = "none",
                         panel.spacing   = unit(0, "lines"),
                         plot.margin     = margin(0, 0, 0, 0)),
               scale_x_continuous(expand = c(0, 0)),
               scale_y_continuous(expand = c(0, 0)))

alltheme   <- theme_bw() + theme(text = element_text(size = 23),
              axis.line       = element_line(color = "black", size = 1),
              axis.ticks      = element_line(color = "black", size = 1),
              panel.border    = element_blank(), panel.grid = element_blank(),
              legend.position = "none")
col_scheme <- make_col_scheme(chars = c("A", "C", "G", "T"),
             cols = c("#5CC93B", "#0D00C4", "#F4B63F", "#BB261A"))
font       <- "helvetica_bold"

### global
letter_order      <- NULL
letter_complement <- NULL

### parse JSON file
parse_json <- function(file, verbose=FALSE) {
  parsed <- rjson::fromJSON(file = file)

  if (verbose) {
    cat("Getting model data\n")
  }

  order <- parsed$modelSettings$letterOrder
  if (is.null(order)) {
    order <- "ACGT"
  }

  complement <- parsed$modelSettings$letterComplement
  if (is.null(complement)) {
    complement <- "A-T,C-G"
  }

  get_complement(complement, order)

  if (!setequal(letter_order, c("A", "C", "G", "T"))) {
    col_scheme <<- "auto"
    font <<- "helvetica_regular"
  }

  if (verbose) {
    cat("Getting coefficients\n")
  }

  coeff <- list()
  for (mode in seq_along(parsed$modelSettings$bindingModes)) {
    coeff[[mode]] <- list()

    coeff[[mode]][["mononucleotide"]] <- mononucleotide_matrix(
      parsed$coefficients$bindingModes[[mode]]$mononucleotide)

    coeff[[mode]][["dinucleotide"]]   <- dinucleotide_matrix(
      parsed$coefficients$bindingModes[[mode]]$dinucleotide)
  }

  return(coeff)
}

### parse command-line data
get_complement <- function(complement, order) {
  letter_order <<- as.vector(strsplit(order, "")[[1]])
  complement <- do.call(rbind, strsplit(strsplit(complement, ",")[[1]], "-"))

  letter_complement <<- vector()
  for (i in seq_along(letter_order)) {
    pos <- which(complement == letter_order[i], arr.ind = TRUE)
    pos[1, 2] <- if (pos[1, 2] == 2) 1 else 2
    letter_complement[i] <<- complement[pos][1]
  }
}

mononucleotide_matrix <- function(vect) {
  if (length(unique(vect)) <= 1) {
    return(NULL)
  }

  vect <- matrix(vect, nrow = length(letter_order))
  rownames(vect) <- letter_order

  return(vect)
}

dinucleotide_matrix <- function(vect) {
  if (length(vect) == 0) {
    return(NULL)
  }

  size <- length(vect[[1]]) / (length(letter_order)^2) + 1

  vect_list <- list()
  for (distance in seq_along(vect)) {
    vect_list[[distance]] <- array(vect[[distance]],
      c(length(letter_order), length(letter_order), size - distance),
      list(letter_order, letter_order, 1:(size - distance)))
  }

  return(vect_list)
}

### logo generation
a <- "
get_info_matrix <- function(pwm) {
  pwm <- exp(pwm)
  bit <- log2(nrow(pwm))
  pwm <- apply(pwm, 2, function(column) column / sum(column))
  out <- apply(pwm, 2,
    function(column) column * sum(bit + log2(column[column > 0])))
  return(pwm)
}
"

mononucleotide_logo <- function(pwm, type="energy", axes=TRUE, reverse=FALSE,
  labels=TRUE) {
  # Centering
  pwm <- apply(pwm, 2, function(column) column - mean(column))

  if (type == "info") {
    pwm <- exp(pwm)
  }

  if (type == "prob") {
    pwm <- apply(exp(pwm), 2, function(column) column / sum(column))
  }

  if (reverse) {
    pwm <- pwm[, rev(seq_len(ncol(pwm)))]

    if (is.null(dim(pwm))) {
      pwm <- matrix(pwm, length(pwm), 1)
    }

    rownames(pwm) <- letter_complement
  }

  if (type == "energy") {
    plot <- ggseqlogo(pwm, method = "c", font = font, col_scheme = col_scheme) +
      list(labs(x = NULL, y = NULL),
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
        alpha = 0.5, fill = "white"),
      scale_y_continuous(breaks = pretty_breaks()),
      geom_hline(yintercept = 0), alltheme)

    if (labels) {
      plot <- plot + ylab(expression(paste(Delta, Delta, "G/RT")))
    }
  }

  if (type == "info") {
    plot <- ggseqlogo(pwm, method = "b", font = font, col_scheme = col_scheme) +
      list(labs(x = NULL, y = NULL),
      scale_y_continuous(breaks = pretty_breaks()), alltheme)

    if (labels) {
      plot <- plot + ylab("Bits")
    }
  }

  if (type == "prob") {
    plot <- ggseqlogo(pwm, method = "c", font = font, col_scheme = col_scheme) +
      list(labs(x = NULL, y = NULL),
      scale_y_continuous(breaks = c(0.0, 0.5, 1.0), limits = c(0.0, 1.0)),
      alltheme)

    if (labels) {
      plot <- plot + ylab("Probability")
    }
  }

  suppressMessages(
    if (!axes) {
      plot <- plot + no_axes
    }
    else {
      plot <- plot + scale_x_continuous(expand = c(0.01, 0),
        breaks = 1:dim(pwm)[2])
    }
  )

  return(plot)
}

dinucleotide_logo <- function(dinucleotide, mononucleotide=NULL, axes=TRUE,
  reverse=FALSE) {
  size   <- dim(dinucleotide[[1]])[3] + 1
  width  <- dim(dinucleotide[[1]])[1]
  values <- matrix(NA, size * width, size * width,
    dimnames = list(1:(size * width), 1:(size * width)))

  if (!is.null(mononucleotide)) {
    for (x in 0:(size - 1)) {
      for (a in 1:width) {
        values[width * x + a, width * x + a] <- mononucleotide[a, x + 1]
      }
    }
  }

  for (dist in seq_along(dinucleotide)) {
    for (x in 0:(size - dist - 1)) {
      for (a in 1:width) {
        for (b in 1:width) {
          values[width * x + b, width * (x + dist) + a] <-
            dinucleotide[[dist]][a, b, x + 1]
          values[width * (x + dist) + a, width * x + b] <-
            dinucleotide[[dist]][a, b, x + 1]
        }
      }
    }
  }

  if (reverse) {
    rownames(values) <- rev(rownames(values))
    colnames(values) <- rev(colnames(values))
  }

  plot <- ggplot(type.convert(as.data.frame.table(values)),
    aes(x = Var1, y = Var2)) +
    list(coord_fixed(ratio = 1), alltheme,
      scale_fill_gradient2(low = "blue", mid = "white", high = "red"))

  # add lines
  if (axes) {
    line1 <- data.frame(x = c(0, size) * width + 0.5,
      y = rep(0:size, each = 2) * width + 0.5)

    line2 <- data.frame(x = rep(0:size, each = 2) * width + 0.5,
      y = c(0, size) * width + 0.5)

    plot  <- plot + list(geom_tile(color = "black", aes(fill = Freq)),
      geom_line(data = line1, aes(x = x, y = y, group = y), lwd = 0.2),
      geom_line(data = line2, aes(x = x, y = y, group = x), lwd = 0.2))
  } else {
    plot  <- plot + list(geom_tile(aes(fill = Freq)))
  }

  # add axis
  if (axes) {
    plot <- suppressMessages(plot +
      scale_x_continuous(expand = c(0.01, 0), breaks =
        seq(floor(width / 2) + 0.5, size * width, width), labels = 1:size) +
      scale_y_reverse(expand = c(0.01, 0), breaks =
        seq(floor(width / 2) + 0.5, size * width, width), labels = 1:size) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()))
  } else {
    plot <- suppressMessages(plot + no_axes + scale_y_reverse(expand = c(0, 0)))
  }

  #plot(plot, vp = grid::viewport(width=0.7, height=0.7, angle=45))
  return(plot)
}

### main function
main <- function(model, out_file, type="energy", mode=1, feat="mono", axes=TRUE,
  rev=FALSE, verbose=FALSE) {
  model  <- parse_json(model)
  if (! (type %in% c("energy", "prob", "info"))) {
    cat("Plot type must be either 'energy', 'prob', or 'info'")
  }

  if (! (feat %in% c("mono", "di", "both"))) {
    cat("Plot features must be either 'mono', 'di', or 'both'")
  }

  height <- 4
  mononucleotide <- model[[mode + 1]][["mononucleotide"]]
  dinucleotide   <- model[[mode + 1]][["dinucleotide"]]

  rev_str  <- if (rev) "rev" else "for"
  axes_str <- if (axes) ""   else "_na"

  if ((feat == "mono" || feat == "both") && !is.null(mononucleotide)) {
    width <- dim(mononucleotide)[2]
    name  <- paste(out_file, "_bm", mode, "_mo_", rev_str, "_", type, axes_str,
      ".png", sep = "")

    if (verbose) {
      cat(paste("Creating", name, "\n"))
    }

    monoplot <- mononucleotide_logo(mononucleotide,
      type = type, reverse = rev, axes = axes)
    ggsave(name, plot = monoplot, width = width, height = height,
      type = "cairo", limitsize = FALSE)
  }

  if ((feat == "di" || feat == "both") && !is.null(dinucleotide)) {
    width <- length(dinucleotide[[1]]) / (length(letter_order)^2) + 1
    name  <- paste(out_file, "_bm", mode, "_di_", rev_str, axes_str, ".png",
      sep = "")

    if (verbose) {
      cat(paste("Creating", name, "\n"))
    }

    diplot <- dinucleotide_logo(dinucleotide, mononucleotide,
      reverse = rev, axes = axes)
    ggsave(name, plot = diplot, width = width, height = width, type = "cairo",
      limitsize = FALSE)
  }

  if (feat == "both" && !is.null(mononucleotide) && !is.null(dinucleotide)) {
    name <- paste(out_file, "_bm", mode, "_bo_", rev_str, "_", type, axes_str,
      ".png", sep = "")

    if (verbose) {
      cat(paste("Creating", name, "\n"))
    }

    both <- plot_grid(monoplot, diplot, ncol = 1, align = "v", axis = "trbl",
      rel_heights = c(height, width))
    ggsave(name, plot = both, width = width, height = height + width,
      type = "cairo", limitsize = FALSE)
  }
}



