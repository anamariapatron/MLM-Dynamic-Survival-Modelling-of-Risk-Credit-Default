
# ===============================
#  testing: we take the mean and approximate the risk probability  
# ===============================

matrices <- list(x1_matrix, x2_matrix, x3_matrix, x4_matrix)

get_percentiles <- function(mat, probs = c(0.5)) {
  apply(mat, 2, function(col) quantile(col, probs = probs))
}

percentiles_list <- lapply(matrices, get_percentiles)
percentiles_matrix <- do.call(rbind, percentiles_list)

betas <- matrix(c(
  50,  0.9,  0.8, -0.2,
  51, 0.92, 0.81, -0.21,
  49, 0.88, 0.79, -0.19,
  52, 0.95, 0.82, -0.22,
  48, 0.89, 0.78, -0.18,
  5,  0.9,  0.8, -0.2,
  51, 0.91, 0.81, -0.21
), nrow = 7, byrow = TRUE)
# Interest rate (+), # Inflation (+),# LTI (+), # Age (-)

dot_products <- sapply(1:ncol(percentiles_matrix), function(i) {
  sum( betas[i, ]*percentiles_matrix[, i] )
})

dot_products <- sapply(1:ncol(percentiles_matrix), function(i) {
  sum(betas[i, ] * percentiles_matrix[, i])
})

dot_products <- matrix(dot_products, ncol = 1)
dot_products


#median values
dot_products_exp <- exp(dot_products)
dot_products_exp


weibull_vals <- mapply(function(shape, scale) {
  rweibull(1, shape = shape, scale = scale)
}, shape = thetas[,1], scale = thetas[,2])

weibull_vals <- matrix(weibull_vals, ncol = 1)
weibull_vals


result <- weibull_vals * dot_products_exp
result

#Conclusion: the risk of the event increases, landmarking is not that accurate but aft is even worse
# [,1]
# [1,] 0.009825581
# [2,] 0.002075948
# [3,] 0.033560692
# [4,] 0.647176944
# [5,] 0.098796546
# [6,] 0.381715775
# [7,] 0.284256476

