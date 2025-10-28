# ballpark Pcod total biomass in 2238:
b <- 322121.7577 #mt
# go to mgn
b <- b * 1e+9 / 20 / 5.7
biom <- seq(from = 0, to = b*10, length.out = 100)
k <- 0.3
sp <- biom * k
alpha <- 4100962763
beta_original <- 1.65e+11

# Original curve
r = (sp * alpha) / (biom + beta_original)
plot(r ~ biom, type = "l", col = "black", lwd = 2)
abline(v = b)
abline(v = b/3, col = "red")

# Add curves for beta multiplied by 2,3,4,5,6,7,8,9,10
colors <- rainbow(19)
for(i in 2:30) {
  beta_new <- beta_original * i
  r_new = (sp * alpha) / (biom + beta_new)
  lines(r_new ~ biom, col = colors[i-1], lwd = 2)
}

# Add legend
legend("topright", 
       legend = c("Original", paste("Beta x", 2:20)),
       col = c("black", colors),
       lwd = 2)

