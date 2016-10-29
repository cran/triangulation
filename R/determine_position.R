#'Determine position of observer
#'
#'Determine the position of an observer based on angles between three known points as seen by the observer.
#'At least two angles must be provided - preferably \code{observer_angle_AB} and \code{observer_angle_AC}
#'(since this combination allows for solutions outside the triangle formed by the points \code{A}, \code{B} and \code{C})
#'
#' @param A A point defined by a vector containing an x- and an y-coordinate
#' @param B A point defined by a vector containing an x- and an y-coordinate
#' @param C A point defined by a vector containing an x- and an y-coordinate
#' @param observer_angle_AB An angle (numeric) expressed in radians (or alternatively the symbol \code{NA})
#' @param observer_angle_AC An angle (numeric) expressed in radians (or alternatively the symbol \code{NA})
#' @param observer_angle_BC An angle (numeric) expressed in radians (or alternatively the symbol \code{NA})
#' @param output_plot Boolean variable indicating whether a plot should be created
#' @param lines_in_plot Boolean variable indicating whether lines should be drawn in the plot
#' @param coordinates_in_plot Boolean variable indicating whether the coordinates should be printet in the plot
#' @param decimals_in_plot Integer indicating the number of decimals used
#'
#' @return Coordinates indicating the observers position. Note that several solutions might exist.
#' @export determine_position
#' @importFrom graphics plot plot.new points segments text
#'
#' @examples
#' determine_position(A = c(0, 0), B = c(10, 0), C = c(5, 5 * 3^0.5), observer_angle_AB = pi * 2/3,
#' observer_angle_AC = pi * 1/2)
#'
#' determine_position(A = c(0, 0), B = c(10, 0), C = c(5, 5), observer_angle_AB = pi * 5/6,
#' observer_angle_AC = pi * 1/2, observer_angle_BC = NA, lines_in_plot = FALSE)
#'
#' determine_position(A = c(0, 0), B = c(10, 0), C = c(5, 5), observer_angle_AB = pi * 5/6,
#' observer_angle_AC = pi * 1/2, observer_angle_BC = pi * 2/3, lines_in_plot = FALSE)

determine_position <- function(A, B, C, observer_angle_AB, observer_angle_AC, observer_angle_BC = NA, output_plot = TRUE, lines_in_plot = TRUE,
    coordinates_in_plot = TRUE, decimals_in_plot = 2) {

    if (length(A) != 2 | length(B) != 2 | length(C) != 2 | is.numeric(A) == FALSE | is.numeric(B) == FALSE | is.numeric(C) == FALSE) {
        stop("An x and y coordinate must be provided for the points A, B and C")
    }

    if (is.numeric(observer_angle_AB) + is.numeric(observer_angle_AC) + is.numeric(observer_angle_BC) <= 1) {
        stop("At least two angles must be provided as numeric - preferably 'observer_angel_AB' and 'observer_angel_AC'")
    }

    # Define known points A, B, C
    x_coords <- c(A[1], B[1], C[1])
    y_coords <- c(A[2], B[2], C[2])

    # Calculate 'observer_angle_AB' or 'observer_angle_AB' if NA and possible - assuming observer is located inside triangle
    if (is.numeric(observer_angle_BC)) {
        if (is.na(observer_angle_AB)) {
            observer_angle_AB <- 2 * pi - observer_angle_AC - observer_angle_BC
            warning("Assuming observer is inside triangle")
        }

        if (is.na(observer_angle_AC)) {
            observer_angle_AC <- 2 * pi - observer_angle_AB - observer_angle_BC
            warning("Assuming observer is inside triangle")
        }
    }

    # Calculate length two sides of triangle
    length_AB <- sqrt((x_coords[2] - x_coords[1])^2 + (y_coords[2] - y_coords[1])^2)
    length_AC <- sqrt((x_coords[3] - x_coords[1])^2 + (y_coords[3] - y_coords[1])^2)

    # Shift such that A i located at (0;0)
    coords_A <- matrix(c(x_coords - x_coords[1], y_coords - y_coords[1]), nrow = 2, byrow = TRUE)

    # Create matrices for transformation and change basis
    P_AB <- cbind(coords_A[, 2], matrix(c(0, 1, -1, 0), ncol = 2) %*% coords_A[, 2, drop = FALSE])/length_AB
    P_AC <- cbind(coords_A[, 3], matrix(c(0, 1, -1, 0), ncol = 2) %*% coords_A[, 3, drop = FALSE])/length_AC

    coords_AB <- solve(P_AB) %*% coords_A
    coords_AC <- solve(P_AC) %*% coords_A

    # Calculate points in new orthogonal basis AB
    peri_coords_AB <- c(length_AB/2, length_AB/2 * tan(pi/2 - observer_angle_AB/2))

    if (peri_coords_AB[2] == 0) {
        center_coords_AB <- peri_coords_AB
    } else {
        center_coords_AB <- c(peri_coords_AB[1], peri_coords_AB[2]/2 - peri_coords_AB[1]^2/peri_coords_AB[2]/2)
    }
    center_coords_AB <- cbind(center_coords_AB, c(center_coords_AB[1], -center_coords_AB[2]))

    # Calculate points in new orthogonal basis AC
    peri_coords_AC <- c(length_AC/2, length_AC/2 * tan(pi/2 - observer_angle_AC/2))

    if (peri_coords_AC[2] == 0) {
        center_coords_AC <- peri_coords_AC
    } else {
        center_coords_AC <- c(peri_coords_AC[1], peri_coords_AC[2]/2 - peri_coords_AC[1]^2/peri_coords_AC[2]/2)
    }
    center_coords_AC <- cbind(center_coords_AC, c(center_coords_AC[1], -center_coords_AC[2]))

    # Transform points back to orthonormal basis
    center_coords_A <- cbind(P_AB %*% center_coords_AB, P_AC %*% center_coords_AC)
    center_coords <- center_coords_A + c(x_coords[1], y_coords[1])

    # Calculate slope of line though ('positive' and 'negative') centers of circles
    center_slope_pp <- center_coords_A[, 1] - center_coords_A[, 3]
    center_slope_pn <- center_coords_A[, 1] - center_coords_A[, 4]
    center_slope_np <- center_coords_A[, 2] - center_coords_A[, 3]
    center_slope_nn <- center_coords_A[, 2] - center_coords_A[, 4]

    # Create transition matrices and calculate center coords i new basis
    P_pp <- cbind(center_slope_pp, matrix(c(0, 1, -1, 0), ncol = 2) %*% center_slope_pp)/sum(center_slope_pp^2)^0.5
    P_pn <- cbind(center_slope_pn, matrix(c(0, 1, -1, 0), ncol = 2) %*% center_slope_pn)/sum(center_slope_pn^2)^0.5
    P_np <- cbind(center_slope_np, matrix(c(0, 1, -1, 0), ncol = 2) %*% center_slope_np)/sum(center_slope_np^2)^0.5
    P_nn <- cbind(center_slope_nn, matrix(c(0, 1, -1, 0), ncol = 2) %*% center_slope_nn)/sum(center_slope_nn^2)^0.5

    centers_pp <- solve(P_pp) %*% center_coords_A[, c(1, 3)]
    centers_pn <- solve(P_pn) %*% center_coords_A[, c(1, 4)]
    centers_np <- solve(P_np) %*% center_coords_A[, c(2, 3)]
    centers_nn <- solve(P_nn) %*% center_coords_A[, c(2, 4)]

    # Calculate circle intercepts in orthonormal basis
    intercept_pp_orthonorm <- P_pp %*% c(0, 2 * centers_pp[2, 2]) + c(x_coords[1], y_coords[1])
    intercept_pn_orthonorm <- P_pn %*% c(0, 2 * centers_pn[2, 2]) + c(x_coords[1], y_coords[1])
    intercept_np_orthonorm <- P_np %*% c(0, 2 * centers_np[2, 2]) + c(x_coords[1], y_coords[1])
    intercept_nn_orthonorm <- P_nn %*% c(0, 2 * centers_nn[2, 2]) + c(x_coords[1], y_coords[1])

    # Determine correct intercept
    test_intercepts <- function(intercept) {
        lines_from_intercept <- matrix(c(x_coords, y_coords), nrow = 2, byrow = TRUE) - matrix(intercept, nrow = 2, ncol = 3)
        dist_to_points <- colSums(lines_from_intercept^2)^0.5

        if (dist_to_points[1] > 10^(-decimals_in_plot) & dist_to_points[2] > 10^(-decimals_in_plot) & dist_to_points[3] > 10^(-decimals_in_plot)) {
            if (abs(sum(lines_from_intercept[, 1] * lines_from_intercept[, 2])/(dist_to_points[1] * dist_to_points[2]) - cos(observer_angle_AB)) <
                1e-06) {
                if (abs(sum(lines_from_intercept[, 1] * lines_from_intercept[, 3])/(dist_to_points[1] * dist_to_points[3]) - cos(observer_angle_AC)) <
                  1e-06) {
                  if (is.na(observer_angle_BC) | abs(sum(lines_from_intercept[, 2] * lines_from_intercept[, 3])/(dist_to_points[2] *
                    dist_to_points[3]) - cos(observer_angle_BC)) < 1e-06) {
                    return(intercept)
                  }
                }
            }
        }
    }

    correct_intercept <- c(test_intercepts(intercept_pp_orthonorm), test_intercepts(intercept_pn_orthonorm), test_intercepts(intercept_np_orthonorm),
        test_intercepts(intercept_nn_orthonorm))

    # Test if solution was found
    if (is.null(correct_intercept)) {
        if (output_plot == TRUE) {
            spread <- max(diff(range(x_coords)), diff(range(y_coords)))/16
            plot(x = x_coords, y = y_coords, ylim = c(min(y_coords) - spread, max(y_coords) + spread), xlab = "x-coordinates",
                ylab = "y-coordinates", asp = 1)
            text(x = x_coords, y = y_coords + spread * 0.75, labels = c("A", "B", "C"), cex = 0.9)
        }
        warning("No solution was found")
        return(c(NULL))
    } else {
        # Tidying up output
        observer_position <- matrix(round(correct_intercept, decimals_in_plot + 6), nrow = 2)
        observer_position <- matrix(unlist(unique(split(observer_position, rep(1:ncol(observer_position), each = nrow(observer_position))))),
            nrow = 2)
        row.names(observer_position) <- c("x-coordinate", "y-coordinate")

        if (output_plot == TRUE) {
            # Create plot
            spread <- max(diff(range(x_coords, observer_position[1, ])), diff(range(y_coords, observer_position[2, ])))/16
            plot(x_coords, y_coords, asp = 1, pch = 1, cex = 1.2, main = "Position of observer", xlim = c(min(x_coords, observer_position[1,
                ]) - spread, max(x_coords, observer_position[1, ]) + spread), ylim = c(min(y_coords, observer_position[2, ]) -
                spread, max(y_coords, observer_position[2, ]) + spread), xlab = "x-coordinates", ylab = "y-coordinates")
            text(x = x_coords, y = y_coords + spread * 0.75, labels = c("A", "B", "C"), cex = 0.9)
            points(x = observer_position[1, ], y = observer_position[2, ], pch = 16, cex = 1.2)

            if (coordinates_in_plot == TRUE) {
                text(x = observer_position[1, ], y = observer_position[2, ] + spread * 0.75, paste0("(", round(observer_position[1,
                  ], decimals_in_plot), " ; ", round(observer_position[2, ], decimals_in_plot), ")"), cex = 0.9)
            }

            if (lines_in_plot == TRUE) {
                segments(x0 = rep(observer_position[1, ], each = 3), y0 = rep(observer_position[2, ], each = 3), x1 = x_coords,
                  y1 = y_coords, lty = "dotted")
            }
        }

        # Return solution
        return(observer_position)
    }
}
