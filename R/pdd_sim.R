pdd_update_lamu <- function(lamu, Phi, K, model) {
    if(model == "a") {#linear PD dependence in speciation rate
        newla <- max(0, lamu[1, 1] * (1 - Phi / K))
        newmu <- lamu[1, 2]
    }

    return(c(newla, newmu))
}

pdd_sum_rates <- function(lamu, N, i) {
    return(lamu[i, 1] * N + lamu[i, 2] * N)
}

pdd_sample_event <- function(lamu, N, i) {
    events = c("spec", "ext", "fake_spec", "fake_ext")

    if((lamu[i - 1, 1] - lamu[i, 1]) >= 0) {
        rspec <- lamu[i, 1]
        rfake_spec <- lamu[i - 1, 1] - lamu[i, 1]
    } else {
        
    }
    if((lamu[i - 1, 2] - lamu[i, 2]) >= 0) {
        rext <- lamu[i, 2]
        rfake_ext <- lamu[i - 1, 2] - lamu[i, 2]
    } else {

    }

    rates <- c(rspec, rext, rfake_spec, rfake_ext)

    return(DDD::sample2(events, 1, prob = rates))
}

pdd_sim<- function (pars, age, model = "a", metric = "pd") {
    if(pars[1] < pars[2]) {stop('the function is designed for lambda_0 > mu_0')}

    if(pars[2] < 0) {stop('per species rates should be positive')}

    if(pars[3] < 0) {stop('clade level carrying capacity should be positive')}

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
        K <- pars[3]
        linlist <- c(-1, 2)
        newL <- 2
        lamu <- matrix(c(pars[1], pars[2]), ncol = 2)
        Phi[i] <- 0

        t[i + 1] <- t[i] + stats::rexp(1, pdd_sum_rates(lamu, N, i))

        # main simulation circle
        while (t[i + 1] <= age) {
            i <- i + 1
            ranL <- sample2(linlist, 1)
            Phi[i] <- L2Phi(L, t[i], metric)
            lamu <- rbind(lamu, pdd_update_lamu(lamu, Phi[i], K, model))
            event <- pdd_sample_event(lamu, N, i)
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
                    lamu[i, ] <- pdd_update_lamu(lamu, Phi[i], K, model)
                }
            } else if (event == "fake_spec" | event == "fake_ext") {
                N[i] <- N[i - 1]
            }
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
    lamuphis <- data.frame("time" = t[-i], "lambda" = lamu[, 1], "mu" = lamu[, 2], "Phi" = Phi, "N" = N)
    out <- list(tes = tes, tas = tas, L = L, brts = brts, lamuphis = lamuphis)

    return(out)
}