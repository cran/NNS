#' NNS Partition Map
#'
#' Creates partitions based on partial moment quadrant centroids, iteratively assigning identifications to observations based on those quadrants (unsupervised partitional and hierarchial clustering method).  Basis for correlation, dependence \link{NNS.dep}, regression \link{NNS.reg} routines.
#'
#' @param x a numeric vector.
#' @param y a numeric vector with compatible dimensions to \code{x}.
#' @param Voronoi logical; \code{FALSE} (default) Displays a Voronoi type diagram using partial moment quadrants.
#' @param type \code{NULL} (default) Controls the partitioning basis.  Set to \code{(type = "XONLY")} for X-axis based partitioning.  Defaults to \code{NULL} for both X and Y-axis partitioning.
#' @param order integer; Number of partial moment quadrants to be generated.  \code{(order = "max")} will institute a perfect fit.
#' @param obs.req integer; (8 default) Required observations per cluster where quadrants will not be further partitioned if observations are not greater than the entered value.  Reduces minimum number of necessary observations in a quadrant to 1 when \code{(obs.req = 1)}.
#' @param min.obs.stop logical; \code{TRUE} (default) Stopping condition where quadrants will not be further partitioned if a single cluster contains less than the entered value of \code{obs.req}.
#' @param noise.reduction the method of determining regression points options for the dependent variable \code{y}: ("mean", "median", "mode", "off"); \code{(noise.reduction = "mean")} uses means for partitions.  \code{(noise.reduction = "median")} uses medians instead of means for partitions, while \code{(noise.reduction = "mode")} uses modes instead of means for partitions.  Defaults to \code{(noise.reduction = "off")} where an overall central tendency measure is used, which is the default for the independent variable \code{x}.
#' @return Returns:
#'  \itemize{
#'   \item{\code{"dt"}} a \code{data.table} of \code{x} and \code{y} observations with their partition assignment \code{"quadrant"} in the 3rd column and their prior partition assignment \code{"prior.quadrant"} in the 4th column.
#'   \item{\code{"regression.points"}} the \code{data.table} of regression points for that given \code{(order = ...)}.
#'   \item{\code{"order"}}  the \code{order} of the final partition given \code{"min.obs.stop"} stopping condition.
#'   }
#'
#' @note \code{min.obs.stop = FALSE} will not generate regression points due to unequal partitioning of quadrants from individual cluster observations.
#'
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments" (ISBN: 1490523995)
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rnorm(100) ; y <- rnorm(100)
#' NNS.part(x, y)
#'
#' ## Data.table of observations and partitions
#' NNS.part(x, y, order = 1)$dt
#'
#' ## Regression points
#' NNS.part(x, y, order = 1)$regression.points
#'
#' ## Voronoi style plot
#' NNS.part(x, y, Voronoi = TRUE)
#'
#' ## Examine final counts by quadrant
#' DT <- NNS.part(x, y)$dt
#' DT[ , counts := .N, by = quadrant]
#' DT
#' }
#' @export


NNS.part = function(x, y,
                    Voronoi = FALSE,
                    type = NULL,
                    order = NULL,
                    obs.req = 8,
                    min.obs.stop = TRUE,
                    noise.reduction = "off"){

    noise.reduction <- tolower(noise.reduction)
    if (!any(noise.reduction %in% c("mean", "median", "mode",
                                    "off", "mode_class"))) {
        stop("Please ensure noise.reduction is from 'mean', 'median', 'mode' or 'off'")
    }

    if(any(class(x)%in%c("tbl","data.table"))) x <- as.vector(unlist(x))
    if(any(class(y)%in%c("tbl","data.table"))) y <- as.vector(unlist(y))

    if (is.null(obs.req)) obs.req <- 8

    if (!is.null(order) && order == 0) order <- 1


    if (Voronoi) {
        x.label <- deparse(substitute(x))
        y.label <- deparse(substitute(y))
    }

    x <- as.numeric(x)
    y <- as.numeric(y)

    PART <- data.table::data.table(x, y, quadrant = "q", prior.quadrant = "pq")[, `:=`(counts, .N), by = "quadrant"][, `:=`(old.counts, .N), by = "prior.quadrant"]

    if(Voronoi) plot(x, y, col = "steelblue", cex.lab = 1.5, xlab = x.label, ylab = y.label)

    if (length(x) <= 8) {
        if(is.null(order)){
            order <- 1
            hard.stop <- max(ceiling(log(length(x), 2)), 1)
        } else {
            obs.req <- 0
            hard.stop <- length(x)
        }
    }

    if(is.null(order)) order <- max(ceiling(log(length(x), 2)), 1)

    if(!is.numeric(order)){
        obs.req <- 0
        hard.stop <- max(ceiling(log(length(x), 2)), 1) + 2
    } else {
        obs.req <- obs.req
        hard.stop <- 2*max(ceiling(log(length(x), 2)), 1) + 2
    }

    if(is.null(type)){
        i <- 0L
        while (i >= 0) {
            if(i == order || i == floor(log(length(x), 2))) break

            PART[counts >= obs.req, `:=`(counts, .N), by = quadrant]
            obs.req.rows <- PART[counts >= obs.req, which = TRUE]
            if(length(obs.req.rows)==0) break
            
            PART[old.counts >= obs.req, `:=`(old.counts, .N), by = prior.quadrant]
            old.obs.req.rows <- PART[old.counts >= obs.req, which = TRUE]
            if(min.obs.stop && obs.req > 0 && (length(obs.req.rows) < length(old.obs.req.rows))) break
            
            l.PART <- max(PART$counts)

            if(noise.reduction == "off") {
                if(Voronoi) {
                    if(l.PART > obs.req) {
                        PART[obs.req.rows, {
                            segments(min(x), gravity(y), max(x), gravity(y),
                                     lty = 3)
                            segments(gravity(x), min(y), gravity(x), max(y),
                                     lty = 3)
                        }, by = quadrant]
                    }
                }
                 RP <- PART[obs.req.rows, mapply(function(f,z) as.list(f(z)), list(x = gravity, y = gravity), .SD), by = quadrant,
                           .SDcols = x:y]


            }

            if(noise.reduction == "mean") {
                if(Voronoi) {
                    if(l.PART > obs.req) {
                        PART[obs.req.rows, {
                            segments(min(x), mean(y), max(x), mean(y),
                                     lty = 3)
                            segments(gravity(x), min(y), gravity(x), max(y),
                                     lty = 3)
                        }, by = quadrant]
                    }
                }
                RP <- PART[obs.req.rows, mapply(function(f,z) as.list(f(z)), list(x = gravity, y = mean), .SD), by = quadrant,
                           .SDcols = x:y]
            }



            if(noise.reduction == "median") {
                if(Voronoi) {
                    if(l.PART > obs.req) {
                        PART[obs.req.rows, {
                            segments(min(x), median(y), max(x), median(y),
                                     lty = 3)
                            segments(gravity(x), min(y), gravity(x), max(y),
                                     lty = 3)
                        }, by = quadrant]
                    }
                }
                RP <- PART[obs.req.rows, mapply(function(f,z) as.list(f(z)), list(x = gravity, y = median), .SD), by = quadrant,
                           .SDcols = x:y]
            }

            if(noise.reduction == "mode") {
                if(Voronoi) {
                    if(l.PART > obs.req) {
                        PART[obs.req.rows, {
                            segments(min(x), mode(y), max(x), mode(y),
                                     lty = 3)
                            segments(gravity(x), min(y), gravity(x), max(y),
                                     lty = 3)
                        }, by = quadrant]
                    }
                }

                RP <- PART[obs.req.rows, mapply(function(f,z) as.list(f(z)), list(x = gravity, y = mode), .SD), by = quadrant,
                           .SDcols = x:y]
            }

            if(noise.reduction == "mode_class") {
                if(Voronoi) {
                    if(l.PART > obs.req) {
                        PART[obs.req.rows, {
                            segments(min(x), mode_class(y), max(x), mode_class(y),
                                     lty = 3)
                            segments(gravity_class(x), min(y), gravity_class(x), max(y),
                                     lty = 3)
                        }, by = quadrant]
                    }
                }
                RP <- PART[obs.req.rows, mapply(function(f,z) as.list(f(z)), list(x = gravity, y = mode_class), .SD), by = quadrant,
                           .SDcols = x:y]
            }

                RP[, `:=`(prior.quadrant, (quadrant))]

                PART[obs.req.rows, `:=`(prior.quadrant, (quadrant))]

                PART[RP, on = .(quadrant), `:=`(q_new, {
                    lox = x.x <= i.x
                    loy = x.y <= i.y
                    1L + lox + loy * 2L
                })]

                PART[obs.req.rows, `:=`(quadrant, paste0(quadrant, q_new))]

            if((min(PART$counts) <= obs.req) && i > 0) break
            if(nrow(PART) > length(x)) break
            i = i + 1L
        }

        if(!exists("RP"))  RP <- PART[, c("quadrant", "x", "y")]
        if(!is.numeric(order) || is.null(dim(RP)))  RP <- PART[, c("quadrant", "x", "y")] else RP[, `:=`(prior.quadrant = NULL)]

        PART[, `:=`(counts = NULL, old.counts = NULL, q_new = NULL)]

        RP <- data.table::setorder(RP[], quadrant)[]

        if (Voronoi) {
            title(main = paste0("NNS Order = ", i), cex.main = 2)
            if(min.obs.stop) points(RP$x, RP$y, pch = 15, lwd = 2, col = "red")
        }

        if(min.obs.stop == FALSE) RP <- NULL
        return(list(order = i, dt = PART[], regression.points = RP))
    }

    if(!is.null(type)) {
        i <- 0L
        while (i >= 0) {
            if(nrow(PART) > length(x)) break
            if(i == order || i == floor(log(length(x), 2))) break

            PART[counts > obs.req/2, `:=`(counts, .N), by = quadrant]
            obs.req.rows <- PART[counts > obs.req/2, which = TRUE]
            if(length(obs.req.rows)==0 && obs.req > 0) break
            
            PART[old.counts > obs.req/2, `:=`(old.counts, .N), by = prior.quadrant]
            old.obs.req.rows <- PART[old.counts > obs.req/2, which = TRUE]
            if(obs.req > 0 && (length(obs.req.rows) < length(old.obs.req.rows))) break

            if(noise.reduction == "off") {
                RP <- PART[obs.req.rows, mapply(function(f,z) as.list(f(z)), list(x = gravity, y = gravity), .SD), by = quadrant,
                           .SDcols = x:y]
            }

            if(noise.reduction == "mean") {
                RP <- PART[obs.req.rows, mapply(function(f,z) as.list(f(z)), list(x = gravity, y = mean), .SD), by = quadrant,
                           .SDcols = x:y]
            }

            if(noise.reduction == "mode") {
                RP <- PART[obs.req.rows, mapply(function(f,z) as.list(f(z)), list(x = gravity, y = mode), .SD), by = quadrant,
                           .SDcols = x:y]
            }

            if(noise.reduction == "mode_class") {
                RP <- PART[obs.req.rows, mapply(function(f,z) as.list(f(z)), list(x = gravity, y = mode_class), .SD), by = quadrant,
                           .SDcols=x:y]
            }

            if(noise.reduction == "median") {
                RP <- PART[obs.req.rows, mapply(function(f,z) as.list(f(z)), list(x = gravity, y = median), .SD), by = quadrant,
                           .SDcols = x:y]
            }

                RP[, `:=`(prior.quadrant, (quadrant))]

                PART[obs.req.rows, `:=`(prior.quadrant, (quadrant))]

                PART[RP, on = .(quadrant), `:=`(q_new, {
                    lox = x.x > i.x
                    1L + lox
                })]

                PART[obs.req.rows, `:=`(quadrant, paste0(quadrant, q_new))]


            if((min(PART$counts) <= obs.req) && i > 0) break
            if(nrow(PART) > length(x)) break
            i <- i + 1L
        }

        if(!exists("RP"))  RP <- PART[, c("quadrant", "x", "y")]
        if(!is.numeric(order) || is.null(dim(RP))) RP <- PART[, c("quadrant", "x", "y")] else RP[, `:=`(prior.quadrant = NULL)]

        PART[, `:=`(counts = NULL, old.counts = NULL, q_new = NULL)]

        RP <- data.table::setorder(RP[], quadrant)[]

        if(is.discrete(x)) RP$x <- ifelse(RP$x%%1 < .5, floor(RP$x), ceiling(RP$x))

        if(Voronoi) {
            abline(v = c(PART[ ,min(x), by=prior.quadrant]$V1,max(x)), lty = 3)
            if(min.obs.stop) points(RP$x, RP$y, pch = 15, lwd = 2, col = "red")
            title(main = paste0("NNS Order = ", i), cex.main = 2)
        }

        if(min.obs.stop == FALSE) RP <- NULL
        return(list(order = i, dt = PART[], regression.points = RP))
    }
}
