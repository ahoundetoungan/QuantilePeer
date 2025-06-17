#' @title The QuantilePeer package
#' @description The \pkg{QuantilePeer} package simulates and estimates peer effect models including the quantile-based specification (Houndetoungan, 2025 \doi{10.3982/ECTA21048}), and the models with Constant Elasticity of Substitution (CES)-based social norm (Boucher et al., 2024 \doi{10.3982/ECTA21048}).
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix Matrix
#' @references Boucher, V., Rendall, M., Ushchev, P., & Zenou, Y. (2024). Toward a general theory of peer effects. Econometrica, 92(2), 543-565, \doi{10.3982/ECTA21048}.
#' @references Houndetoungan, A. (2025). Quantile peer effect models. arXiv preprint arXiv:2405.17290, \doi{10.48550/arXiv.2506.12920}.
#' @references Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4), 361-365, \doi{10.1080/00031305.1996.10473566}.
#' @useDynLib QuantilePeer, .registration = TRUE
"_PACKAGE"
NULL