pdd_update_lamu <- function(lamu, Phi, K, N, betas, model) {
    if (is.null(betas)) {
        if (model == "a") {#linear PD dependence in speciation rate
            newla <- max(0, lamu[1, 1] * (1 - Phi / K))
            newmu <- lamu[1, 2]
        }
    }

    if (is.null(K)) {
        if (model == "b") {
            newla <- max(0, lamu[1, 1] + betas[1] * N + betas[2] * Phi)
            newmu <- lamu[1, 2]
        }
    }

    return(c(newla, newmu))
}

pdd_sum_rates <- function(lamu, N, i) {
    return(lamu[i, 1] * N + lamu[i, 2] * N)
}

pdd_sample_event <- function(lamu, N, betas, age, t, i, model) {
    events <- c("spec", "ext", "fake_spec", "fake_ext")
    la_max <- 0
    mu_max <- 0

    if (is.null(betas)) {
        if (model == "a") {
            la_max <- lamu[i - 1, 1]
            mu_max <- lamu[i - 1, 2]

            rspec <- lamu[i, 1]
            rfake_spec <- la_max - lamu[i, 1]

            rext <- lamu[i, 2]
            rfake_ext <- mu_max - lamu[i, 2]
        } 
    } else {
        if (model == "b") {
            if (betas[2] < 0) {
                la_max <- lamu[i, 1]
                mu_max <- lamu[i, 2]

                rspec <- lamu[i + 1, 1]
                rfake_spec <- la_max - lamu[i + 1, 1]

                rext <- lamu[i + 1, 2]
                rfake_ext <- mu_max - lamu[i + 1, 2]
            } else {
                mu_max <- lamu[i, 2]

                rspec <- lamu[i + 1, 1]
                rfake_spec <- (lamu[i + 1, 1] - lamu[i, 1]) / (t[i] - t[i - 1]) * (age - t[i])
                
                rext <- lamu[i + 1, 2]
                rfake_ext <- mu_max - lamu[i + 1, 2]
            }
        } 
    }

    rates <- c(rspec, rext, rfake_spec, rfake_ext)

    return(DDD::sample2(events, 1, prob = rates))
}

pdd_sim<- function (la, mu, K, beta_N, beta_Phi, age, model = "a", metric = "pd") {
    if (missing(K)) {
        if(missing(beta_N) | missing(beta_Phi)) {
            stop('incomplete parameter list')
        }
    } else {
        if(!missing(beta_N) | !missing(beta_Phi)) {
            stop('unused parameters')
        }
    }

    done <- 0
    while (done == 0) {
        # initialization
        t <- rep(0, 1)
        L <- matrix(0, 2, 4)
        i <- 1
        t[1] <- 0
        N <- 2
        L[1, 1 : 4] <- c(0, 0, -1, -1)
        L[2, 1 : 4] <- c(0, -1, 2, -1)
        Phi <- rep(0, 1) # PD
        linlist <- c(-1, 2)
        newL <- 2
        lamu <- matrix(c(la, mu), ncol = 2)
        Phi[i] <- 0

        if (missing(beta_N) | missing(beta_Phi)) {
            betas = NULL
        } else {
            betas <- c(beta_N, beta_Phi)
        }

        if (missing(K)) {
            K = NULL
            lamu <- rbind(lamu, pdd_update_lamu(lamu, Phi[i], K, N, betas, model))
        }

        t[i + 1] <- t[i] + stats::rexp(1, pdd_sum_rates(lamu, N, i))

        # main simulation circle
        while (t[i + 1] <= age) {
            i <- i + 1
            ranL <- sample2(linlist, 1)
            Phi[i] <- L2Phi(L, t[i], metric)
            lamu <- rbind(lamu, pdd_update_lamu(lamu, Phi[i], K, N, betas, model))
            event <- pdd_sample_event(lamu, N, betas, age, t, i, model)

            if (event == "spec") {
                N[i] <- N[i - 1] + 1
                newL <- newL + 1
                L <- rbind(L, c(t[i], ranL, sign(ranL) * newL, 
                  -1))
                linlist <- c(linlist, sign(ranL) * newL)
            } else if (event == "ext") {
                N[i] <- N[i - 1] - 1
                L[abs(ranL), 4] <- t[i]
                w <- which(linlist == ranL)
                linlist <- linlist[-w]
                linlist <- sort(linlist)

                if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
                    
                } else {
                    Phi[i] <- L2Phi(L, t[i], metric)
                }
            } else if (event == "fake_spec" | event == "fake_ext") {
                N[i] <- N[i - 1]
            }

            lamu[i + 1, ] <- pdd_update_lamu(lamu, Phi[i], K, N, betas, model)

            if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
                t[i + 1] <- Inf
            } else {
                t[i + 1] <- t[i] + stats::rexp(1, pdd_sum_rates(lamu, N, i))
            }
        }

        if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
            done <- 0
        } else {
            done <- 1
        }
    }

    L[, 1] <- age - c(L[, 1])
    notmin1 <- which(L[, 4] != -1)
    L[notmin1, 4] <- age - c(L[notmin1, 4])
    L[which(L[, 4] == age + 1), 4] <- -1
    tes <- L2phylo(L, dropextinct = T)
    tas <- L2phylo(L, dropextinct = F)
    brts <- L2brts(L, dropextinct = T)
    lamuphis <- data.frame("time" = t[-i], "lambda" = lamu[-1, 1], "mu" = lamu[-1, 2], "Phi" = Phi, "N" = N)
    out <- list(tes = tes, tas = tas, L = L, brts = brts, lamuphis = lamuphis)

    return(out)
}