spline_profile_matrix <- function(profile_matrix, time, targetTime, df = NULL, knots = NULL, degree = 3) {
  if(any(is.na(profile_matrix))) {
    stop("Profiles cannot be NA")
  }
  #Create the spline basis
  spline_basis = bs(time, df = df, knots = knots, degree = degree)

  #There is a zero-only last column in the basis and it sometimes breaks things, so we have to remove it
  num_coeff <- dim(spline_basis)[2] - 1
  if(any(spline_basis[,num_coeff + 1] != 0)) {
    stop("Broken assumption of zero column in bs output")
  }
  spline_basis_coeffs = spline_basis[,1:num_coeff]

  #Use a least-squares fit of a B-spline basis
  spline_fit <- lm(t(profile_matrix) ~ 0 + spline_basis_coeffs); #Without intercept

  #Create the same basis but for the target time
  basis_new <-  bs(x = targetTime, degree = degree, knots = attr(spline_basis, "knots"), Boundary.knots = attr(spline_basis, "Boundary.knots"));
  basis_new <- basis_new[,1:num_coeff]

  #The result is the product of the spline coefficients with the spline basis for target time
  splined_result = t(basis_new %*% spline_fit$coefficients);

  #Ensure the profiles are strictly positive
  splined_result[splined_result < 0] = 0;

  rownames(splined_result) = rownames(profile_matrix)
  return(splined_result)
}
