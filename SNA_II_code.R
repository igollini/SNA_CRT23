## CODE FOR SNA_II

###
## Networks as Random Graphs
###
library('statnet')
n <- 9
set.seed(1)
y <- network(n, directed = FALSE)
y[,]

plot(y, 
     vertex.cex = 6, 
     vertex.col = 0, 
     label = 1:n, 
     label.pos = 5,
     label.cex = 2)

## ----echo = FALSE, eval = TRUE, message = FALSE, fig.align="center", fig.height = 4, fig.width = 6, out.width = "100%"----
par(mar = c(4, 4, 0, 0))
deg_dist <- table(degree(y, gmode = "graph"))
plot(deg_dist / sum(deg_dist), 
     ylab = 'Frequency',
     xlab = 'Degree',
     lwd = 4)


n_d <- 12
set.seed(1)
y_d <- network(n_d)
y_d[,]

par(mar = rep(0, 4))
plot(y_d, 
     vertex.cex = 5, vertex.col = 0, 
     arrowhead.cex = 3,
     label = 1:n_d, label.pos = 5,
     label.cex = 2)

ideg_dist <- table(degree(y_d, cmode = "indegree"))
odeg_dist <- table(degree(y_d, cmode = "outdegree"))
par(mfrow = c(1, 2), mar = c(4, 4, 0, 4))
plot(ideg_dist/ sum(ideg_dist), 
     xlab = 'In-degree', ylab = 'Frequency', lwd = 4)
plot(odeg_dist/ sum(odeg_dist),
     xlab = 'Out-degree', ylab = 'Frequency', lwd = 4)

###
# Erdos - Renyi
###

e_y <- summary(y ~ edges)
m_y <- choose(n, 2)

eta_mle <- unname(e_y / m_y)
eta_mle

theta_mle <- log(eta_mle / (1 - eta_mle))
theta_mle

eta <- 0.3 
theta <- log(eta / (1 - eta))
n <- 9
set.seed(1)
y_sim <- simulate(~ edges, 
                  coef = theta,
                  basis = network(n, directed = FALSE))

par(mar = rep(0, 4))
plot(y_sim, 
     vertex.cex = 5, vertex.col = 0, 
     label = 1:n, label.pos = 5,
     label.cex = 2)

###
## Dyadic-independence
data(sampson)
y_d <- samplike

par(mar = rep(0, 4))
set.seed(1)
plot(y_d, 
     vertex.cex = 6, vertex.col = 0, 
     label = y_d %v% 'vertex.names', 
     label.pos = 5)

r_model <- ergm(y_d ~ edges + mutual)
summary(r_model)$coefficients

theta_e <- unname(r_model$coefficients["edges"])
theta_r <- unname(r_model$coefficients["mutual"])
z_ij <- 1 + 2 * exp(theta_e) + exp(theta_e + theta_r)
(exp(theta_e) + exp(2 * theta_e + theta_r)) / z_ij

exp(theta_e + theta_r) / (1 + exp(theta_e + theta_r))

par(mar = rep(0, 4))
set.seed(1)
plot(y_d, vertex.cex = 3, 
     vertex.col = unclass(factor(y_d %v% 'group')) + 1)
legend('topright', pt.bg = 2:4, pt.cex = 1.2, pch = 21, 
       legend = sort(unique(y_d %v% 'group')))


r_model_2 <- ergm(y_d ~ edges + mutual + mutual('group'))
summary(r_model_2)$coefficients

edge_prob_y_ij_1 <- function(y_ji, x_ij, theta){
  
  num <- unname(exp(theta["edges"] + 
                      theta["mutual"] * y_ji + 
                      theta["mutual.group"] * y_ji * x_ij))
  den <- 1 + num
  
    return(num / den)
  
}

edge_prob_y_ij_1(y_ji = 0, x_ij = 0, 
                 theta = r_model_2$coefficients)
edge_prob_y_ij_1(y_ji = 1, x_ij = 0, 
                 theta = r_model_2$coefficients)
edge_prob_y_ij_1(y_ji = 0, x_ij = 0, 
                 theta = r_model_2$coefficients)
edge_prob_y_ij_1(y_ji = 1, x_ij = 1, 
                 theta = r_model_2$coefficients)

dyadic_prob <- function(x_ij, theta){
  theta_e <- unname(theta["edges"])
  theta_r <- unname(theta["mutual"])
  theta_h <- ifelse(x_ij == 1, 
                    theta["mutual.group"], 0)
  z_ij <- 1 + 2 * exp(theta_e) + 
    exp(2 * theta_e + theta_r + theta_h)
  
  p_00 <- 1 / z_ij
  p_10 <- exp(theta_e) / z_ij # = p_01
  p_11 <- exp(2 * theta_e + theta_r + theta_h) / z_ij
  
  return(c(p_00 = p_00, p_10 = p_10, p_11 = p_11))
}

plot(dyadic_prob(x_ij = 1, theta = r_model_2$coefficients), 
     type = 'h', xlab = 'Number of edges in the dyad',
     ylab = 'Probability', main = expression(X[ij] == 1),
     lwd = 4, xaxt="n")
axis(1, at = 1:3, labels = 0:2)

par(mar = c(4.5, 4, 4, 0.5))
plot(dyadic_prob(x_ij = 0, 
                 theta = r_model_2$coefficients), 
     type = 'h', 
     xlab = 'Number of edges in the dyad',
     ylab = 'Probability', 
     main = expression(X[ij] == 0), 
     lwd = 4, xaxt="n")
axis(1, at = 1:3, labels = 0:2)

# Conditional odds ratios of edges

y_d[2, 7] <- 1
s_y_ij_1 <- summary(y_d ~ edges + mutual + mutual('group'))

y_d[2, 7] <- 0
s_y_ij_0 <- summary(y_d ~ edges + mutual + mutual('group'))

delta_y_ij <- s_y_ij_1 - s_y_ij_0
delta_y_ij

theta <- unname(r_model_2$coefficients)
OR_y_ij <- c(exp(t(theta) %*% delta_y_ij)) 
OR_y_ij

###
## Markov dependence assumption

par(mar = rep(0, 4))
n <- 4
set.seed(123)
y <- network(n, directed = FALSE)
plot(y, vertex.cex = 4, 
     vertex.col = 0, 
     label = 1:n, label.pos = 5)

sim <- function(theta_triangle){
  simulate(~ edges + kstar(2:3) + triangle,
           basis = network(30, directed = FALSE), coef = c(-3, 0.5, -0.2, theta_triangle), 
           output = 'stats',
           control = control.simulate.formula(MCMC.burnin = 50000, MCMC.interval = 1000),
           nsim = 100)[,4] |>
    quantile(probs = c(0.025, 0.50, 0.975))
}
set.seed(123)
theta_triangle_values <- seq(0, 1.5, by = 0.05)
simt <- sapply(theta_triangle_values, sim)

par(mar = c(4, 4, 1, 0))
plot(x = theta_triangle_values, 
     y = simt[2, ], 
     main = '',
     type = 'l', 
     xlab = 'Triangle parameter',
     ylab = 'Number of triangles')
lines(x = theta_triangle_values, 
      y = simt[1, ],
      lty = 2)
lines(x = theta_triangle_values, 
      y = simt[3, ],
      lty = 2)

triangle.coef <- 1.1502
n.sim <- 1000

set.seed(1)
y_sims <- simulate(~ edges + kstar(2:3) + triangle,
                   basis = network(30, directed = FALSE),
                   coef = c(-3, 0.5, -0.2, triangle.coef), 
                   output = 'stats',
                   control = control.simulate.formula(MCMC.burnin = 50000, 
                                                      MCMC.interval = 1000),
                   nsim = n.sim)

par(mar = c(4, 4, 1, 0))
hist(y_sims[, 4], 
     breaks = 50, 
     col = 'lightblue', 
     main = "",
     xlab = 'Number of triangles')

data(zach)
y <- zach
cp <- hcl.colors(5, 'Tropic') # palette
f_id <- y %v% 'faction.id' + 3
ids <- unique(f_id)
factions <- unique(y %v% 'faction')
par(mar = rep(0, 4)) 
set.seed(1)
plot(y, label = 1:(dim(y[, ])[1]),
     label.cex = 1, label.pos = 5,
     vertex.cex = 4, edge.col = 'grey',
     vertex.col = cp[f_id])
legend('topright', pt.cex = 1.2, pch = 21, 
       pt.bg  = cp[sort(ids)], 
       legend = factions[order(ids)], 
       title  = 'Faction')

par(mar = rep(4, 4))
plot(summary(y ~ esp(0:(n - 1))), 
     type = 'h', 
     ylab = 'Number of edges',
     xlab = expression(esp[k]),
     main = '',
     lwd = 4, 
     xaxt = 'n')
axis(1, at = 1:n, labels = 0:(n - 1))

par(mar = rep(4, 4))
plot(summary(y ~ esp(0:(n - 1))), 
     type = 'h', 
     ylab = 'Number of edges',
     xlab = expression(esp[k]),
     main = '',
     lwd = 4, 
     xaxt = 'n')

axis(1, at = 1:n, labels = 0:(n - 1))

formula <- y ~ edges + nodematch('faction') + gwesp(decay = 0.7, fixed = TRUE)
set.seed(1)
gwm_model <- ergm(formula)
summary(gwm_model)$coefficients

theta <- unname(gwm_model$coefficients)

y[33, 34] <- 1; s_y_ij_1 <- summary(formula)
y[33, 34] <- 0; s_y_ij_0 <- summary(formula)

delta_y_ij <- s_y_ij_1 - s_y_ij_0
delta_y_ij
OR_y_ij <- c(exp(t(theta) %*% delta_y_ij)) 
OR_y_ij

y[10, 17] <- 1
s_y_ij_1 <- summary(formula)

y[10, 17] <- 0
s_y_ij_0 <- summary(formula)

delta_y_ij <- s_y_ij_1 - s_y_ij_0
delta_y_ij
OR_y_ij <- c(exp(t(theta) %*% delta_y_ij)) 
OR_y_ij

y[1, 34] <- 1
s_y_ij_1 <- summary(formula)

y[1, 34] <- 0
s_y_ij_0 <- summary(formula)

delta_y_ij <- s_y_ij_1 - s_y_ij_0
delta_y_ij
OR_y_ij <- c(exp(t(theta) %*% delta_y_ij)) 
OR_y_ij

###
## Estimation and Goodness of fit

par(mar = rep(0, 4))
n <- 20
set.seed(1)
y_sim <- simulate(~ edges + gwesp(0.7, fixed = TRUE), 
                  coef = c(-2, 0.5),
                  basis = network(n, directed = FALSE),
                  control = control.simulate(MCMC.interval = 5000))
plot(y_sim, vertex.cex = 2)

data(zach)
y <- zach
formula <- y ~ edges +  
    nodematch('faction') +
    gwesp(decay = 0.7, fixed = TRUE)

set.seed(1)
gwm_model <- ergm(formula)
gwm_gof <- gof(gwm_model)
par(mfrow = c(2, 2)) 
plot(gwm_gof)

### 
## Latent Space Models

data(zach)
y <- zach
library(latentnet)
system.time(lat.d2 <- 
  ergmm(y ~ euclidean(d = 2)))

plot(lat.d2)

gf.lat.d2 <- gof(lat.d2, 
                 GOF = ~ degree + esp + dsp + triadcensus + distance)
par(mfrow = c(2,3)); plot(gf.lat.d2)


# devtools::install_github("igollini/lvm4net")
library(lvm4net)
system.time(
  lvm.d2 <- lsm(y[,], D = 2, 
                nstart = 5))

plot(lvm.d2, y[,], 
     drawCB = TRUE)

lsmGof <- goflsm(lvm.d2, Y = y[,])

