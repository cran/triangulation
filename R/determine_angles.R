#' Determine angles as seen by observer
#'
#'Determine the angles (between three known points) as seen by an observer with a known position.
#'
#' @param A A point defined by a vector containing an x- and an y-coordinate
#' @param B A point defined by a vector containing an x- and an y-coordinate
#' @param C A point defined by a vector containing an x- and an y-coordinate
#' @param observer_position A vector containing an x- and an y-coordinate
#' @param output_plot Boolean variable indicating whether a plot should be created
#' @param lines_in_plot Boolean variable indicating whether lines should be drawn in the plot
#' @param angles_in_plot Boolean variable indicating whether the angles should be printet in the plot
#' @param decimals_in_plot Integer indicating the number of decimals used
#'
#' @return The angles as seen by the observer expressed in radians.
#' @export determine_angles
#' @importFrom graphics plot plot.new points segments text
#'
#' @examples
#' determine_angles(A = c(0, 0), B = c(10, 0), C = c(5, 5), observer_position=c(4,1))
#'
#' determine_angles(A = c(0, 0), B = c(10, 0), C = c(5, 5), observer_position=c(4,40),
#' angles_in_plot = FALSE)

determine_angles <- function(A, B, C, observer_position = c(0, 0), output_plot = TRUE, lines_in_plot = TRUE, angles_in_plot = TRUE,
    decimals_in_plot = 2) {

    if (length(A) != 2 | length(B) != 2 | length(C) != 2 | length(observer_position) != 2 | is.numeric(A) == FALSE | is.numeric(B) ==
        FALSE | is.numeric(C) == FALSE | is.numeric(observer_position) == FALSE) {
        stop("An x and y coordinate must be provided for the points 'A', 'B', 'C' and 'oberserver_position'")
    }

    # Define known points A, B, C
    x_coords <- c(A[1], B[1], C[1]) - observer_position[1]
    y_coords <- c(A[2], B[2], C[2]) - observer_position[2]

    coords_matrix <- matrix(c(x_coords, y_coords), nrow = 2, byrow = TRUE)

    # Calculate distances from observer to points
    dist_to_points <- colSums(coords_matrix^2)^0.5

    if (dist_to_points[1] < 10^(-decimals_in_plot) | dist_to_points[2] < 10^(-decimals_in_plot) | dist_to_points[3] < 10^(-decimals_in_plot)) {
        stop("Observer can not be located at the points A B or C")
    }

    # Calculate angles as seen by observer
    observer_angle_AB <- acos(sum(coords_matrix[, 1] * coords_matrix[, 2])/(dist_to_points[1] * dist_to_points[2]))
    observer_angle_AC <- acos(sum(coords_matrix[, 1] * coords_matrix[, 3])/(dist_to_points[1] * dist_to_points[3]))
    observer_angle_BC <- acos(sum(coords_matrix[, 2] * coords_matrix[, 3])/(dist_to_points[2] * dist_to_points[3]))

    observer_angles <- c(observer_angle_AB = observer_angle_AB, observer_angle_AC = observer_angle_AC, observer_angle_BC = observer_angle_BC)

    if (output_plot == TRUE) {
        # Create plot
        x_coords <- x_coords + observer_position[1]
        y_coords <- y_coords + observer_position[2]

        spread <- max(diff(range(x_coords, observer_position[1])), diff(range(y_coords, observer_position[2])))/16
        plot(x_coords, y_coords, asp = 1, pch = 1, cex = 1.2, main = "Angles as seen by observer", xlim = c(min(x_coords, observer_position[1]) -
            spread, max(x_coords, observer_position[1]) + spread), ylim = c(min(y_coords, observer_position[2]) - spread, max(y_coords,
            observer_position[2]) + spread), xlab = "x-coordinates", ylab = "y-coordinates")
        text(x = x_coords, y = y_coords + spread * 0.75, labels = c("A", "B", "C"), cex = 0.75)
        points(x = observer_position[1], y = observer_position[2], pch = 16, cex = 1.2)

        if (angles_in_plot == TRUE) {
            text_x <- c(mean(c(A[1], B[1], observer_position[1], observer_position[1])), mean(c(A[1], C[1], observer_position[1],
                observer_position[1])), mean(c(B[1], C[1], observer_position[1], observer_position[1])))
            text_y <- c(mean(c(A[2], B[2], observer_position[2], observer_position[2])), mean(c(A[2], C[2], observer_position[2],
                observer_position[2])), mean(c(B[2], C[2], observer_position[2], observer_position[2])))
            text(x = text_x, y = text_y, labels = paste0("observer_angle\n_", c("AB", "AC", "BC"), " = ", round(observer_angles/pi,
                decimals_in_plot), "pi"), cex = 0.75)
        }

        if (lines_in_plot == TRUE) {
            segments(x0 = rep(observer_position[1], 3), y0 = rep(observer_position[2], 3), x1 = x_coords, y1 = y_coords, lty = "dotted")
        }
    }

    # Return solution
    return(round(observer_angles, decimals_in_plot + 6))
}
