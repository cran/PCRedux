
res <- inder(x = RAS002[, 1], y = RAS002[, 2])

res_cpD <- summary(res)

index <- c(2,1,3)
res_val <- cbind(cpD = res_cpD[index], 
      rfu = sapply(index, function(i) {
            res[which(res[, "x"] == res_cpD[i]), "y"]}
))

res_val[, "rfu"] <- res_val[, "rfu"] / res_val[2, "rfu"]

vec_a <- c(x = res_val[3, 1] - res_val[2, 1], y = res_val[3, 2] - res_val[2, 2])
vec_b <- c(x = res_val[1, 1] - res_val[2, 1], y = res_val[1, 2] - res_val[2, 2])

lampda <- 360 - acos((vec_a %*% vec_b) / (sum(vec_a^2)^0.5 * sum(vec_b^2)^0.5)) / (pi/180) 
