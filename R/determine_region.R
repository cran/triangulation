#'Determine confidence region for position
#'
#'This function is similar to \code{determine_position()}except for the fact that it is assumed that the angles are subject to measurement error.
#'Hence a confidence region (error 'ellipse') is returned instead of an exact position.
#'
#' @param A A point defined by a vector containing an x- and an y-coordinate
#' @param B A point defined by a vector containing an x- and an y-coordinate
#' @param C A point defined by a vector containing an x- and an y-coordinate
#' @param observer_angle_AB An angle (numeric) expressed in radians
#' @param observer_angle_AC An angle (numeric) expressed in radians
#' @param angle_error A numeric indicating the measurement error in radians
#' @param number_of_points A numeric indicating the number of error points tested
#' @param output_plot Boolean variable indicating whether a plot should be created
#' @param lines_in_plot Boolean variable indicating whether lines should be drawn in the plot
#' @param coordinates_in_plot Boolean variable indicating whether the coordinates should be printet in the plot
#' @param decimals_in_plot Integer indicating the number of decimals used
#'
#' @return Coordinates indicating the outer border of the confidence region. Note that several different regions may exist.
#' @export determine_region
#' @importFrom graphics plot plot.new points segments text
#'
#'
#' @examples
#' determine_region(A = c(0, 0), B = c(10, 0), C = c(5, 5 * 3^0.5), observer_angle_AB = pi * 2/3,
#' observer_angle_AC = pi * 1/2)
#'
#' determine_region(A = c(0, 0), B = c(10, 0), C = c(5, 5), observer_angle_AB = pi * 5/6,
#' observer_angle_AC = pi * 1/2, lines_in_plot = FALSE)

determine_region <- function(A, B, C, observer_angle_AB, observer_angle_AC, angle_error = pi/24, number_of_points = 200, output_plot = TRUE,
    lines_in_plot = FALSE, coordinates_in_plot = FALSE, decimals_in_plot = 2) {

    if (length(A) != 2 | length(B) != 2 | length(C) != 2 | is.numeric(A) == FALSE | is.numeric(B) == FALSE | is.numeric(C) == FALSE) {
        stop("An x and y coordinate must be provided for the points A, B and C")
    }

    if (is.numeric(observer_angle_AB) != TRUE | is.numeric(observer_angle_AC) != TRUE | is.numeric(angle_error) != TRUE | is.numeric(number_of_points) !=
        TRUE) {
        stop("'observer_angle_AB', 'observer_angle_AC', 'angle_error' and 'number_of_points' must be provided as numeric")
    }

    # Define uncorrelated errors
    error_1 <- seq(from = -angle_error, to = angle_error, length = round((number_of_points + 2)/2, 0))
    error_2 <- sqrt(angle_error^2 - error_1^2)
    error_2 <- c(error_2[-1], -error_2[-length(error_1)])
    error_1 <- c(error_1[-1], error_1[-length(error_1)])

    # 'loop' though angles using the 'determine_position'-function
    region_coordinates <- suppressWarnings(mapply(determine_position, observer_angle_AB = c(observer_angle_AB + error_1), observer_angle_AC = c(observer_angle_AC +
        error_2), MoreArgs = list(A = A, B = B, C = C, observer_angle_BC = NA, output_plot = FALSE)))

    if (is.null(region_coordinates)) {
        plot.new()
        stop("No solution exists")
    }

    region_coordinates <- matrix(unlist(region_coordinates), nrow = 2)

    row.names(region_coordinates) <- c("x-coordinate", "y-coordinate")

    if (output_plot == TRUE) {
        # Get plot from 'determine_position'-function
        determine_position(A, B, C, observer_angle_AB, observer_angle_AC, NA, TRUE, lines_in_plot, coordinates_in_plot, decimals_in_plot)
        points(region_coordinates[1, ], region_coordinates[2, ], pch = 16, cex = 0.4)
    }

    # Return solution
    return(round(region_coordinates, decimals_in_plot + 6))

}
