rm(list=ls())

# datadir is the user folder for saving data
datadir <- paste0( getwd(),"/" )  # get automatically current work directory 
# for other folder you need to specify it

# ==============================================
# ==============================================
## Setting parameters
B <- 3000; 
NPERM <- 1000;
PARlist <- list(
  n_vec = c( 10, 20, 50, 100, 500 ),
  mu_vec = c(0,.2,.5,.8),
  sigma_vec = c(1,2,3),
  alpha_vec = c(0,2,10)
)
OUTFILE <- "sim09.rda" 

# ==============================================
# UTILITY FUNCTIONS
# ==============================================
#' @name min_dskew_normal
#' @description Calcola il minimo tra due densità Skew-Normal
#' @note La prima densità ha i parametri fissati \code{xi = 0}, 
#' \code{omega = 1}, \code{alpha = 0}, per cui di fatto è normale standard
#' @param x = x vector
#' @param xi = location parameter
#' @param omega = scale parameter
#' @param alpha = skewness parameter
#' @param plot = logical, if \code{TRUE} produce la rappresentazione grafica delle densità 
#' e dell'area di sovrapposizione
#' @param return.all = logical, if \code{TRUE} restituisce il data set completo delle 
#' densità delle due distribuzioni
min_dskew_normal <- function( x = seq( -5, 5, by = .01 ), xi = 0, omega = 1, alpha = 0, 
                              return.all = FALSE ) {
  
  require( sn )
  y1 <- dsn( x, xi = 0, omega = 1, alpha = 0 )
  y2 <- dsn( x, alpha = alpha, xi = xi, omega = omega )
  dy <- ifelse( y1 < y2, y1, y2 )
  gData <- data.frame( x, y1, y2, dy )  
  
  if (return.all) {
    return( list( gData = gData ) )
  } else {
    return( dy )  
  }
}


# ++++++++++++++++++++++++++++
#' @name permTest
#' @description Esegue test di permutazione su overlapping, 
#'differenza tra medie e rapporto tra varianze
#' @param xList = lista di due elementi (\code{x1} e \code{x2} ) 
#' @param B = numero di permutazioni da effettuare
#' @note Il confronto tra le medie è ad una sola coda e 
#'testa l'ipotesi che le medie siano uguali vs l'ipotesi
#'che la seconda sia maggiore della prima (\code{mean(x2) > mean(x1)})
#' @return Restituisce una lista con tre elementi:
#' obs = vettore dei valori osservati di non-sovrapposizione 
#'       \coed{1-eta}, differenza tra le medie (\code{mean(x2)-mean(x1)}), 
#'       rapporto tra le varianze
#' perm = matrice Bx3 con i valori delle stesse statistiche ottenute
#'        via permutazione
#' pval = vettore con i tre p-values   
permTest <- function( xList, B = 1000) {
  
  require(overlapping)
  names(xList) <- c("x1","x2")
  N <- unlist( lapply(xList,length) )
  
  # observed statistics
  zobs <- 1-overlap( xList )$OV
  dobs <- diff( unlist( lapply(xList, mean) ) )
  Fobs <-  with( xList, var.test(x1,x2)$statistic )
  OBS <- c(zobs,dobs,Fobs)
  names(OBS) <- c("z","d","F")
  Mobs <- matrix( OBS, nrow=B+1, ncol=3, byrow = TRUE )
  
  Yperm <- t( sapply(1:B, function(b){
    xperm <- sample( unlist( xList ) )
    xListperm <- list( x1 = xperm[1:N[1]], x2 = xperm[(N[1]+1):(sum(N))] )
    
    zperm <- 1 - overlap( xListperm )$OV
    dperm <- diff( unlist( lapply(xListperm, mean) ) )
    Fperm <-  with( xListperm, var.test(x1,x2)$statistic )
    
    out <- c(zperm,dperm,Fperm)
    names(out) <- c("zperm","dperm","Fperm")
    out
  }) )
  
  Yperm <- rbind( Mobs[1,], Yperm  )
  
  PVAL <- apply( Yperm > Mobs, 2, mean )
  
  L <- list(obs=OBS,perm=Yperm,pval=PVAL)
  
  return(L)
}


# +++++++++++++++++++++++++++
snpar <- function(xi=0,omega=1,alpha=0) {
  delta <- alpha/sqrt(1+alpha^2)
  mu <- xi + omega * delta * sqrt( 2/pi )
  sigma2 <- omega^2 * ( 1 - (2*delta^2)/pi )
  return(list(mu = mu, sigma = sqrt(sigma2)))
}

# +++++++++++++++++++++++++++
sninvpar <- function( mu=0, sigma=1, xi=NULL, omega=NULL, alpha=0 ) {
  
  if (is.null(omega)) {
    delta <- alpha/sqrt(1+alpha^2)
    omega2 <- sigma^2 / ( 1 - (2*delta^2) / pi )
    omega <- sqrt( omega2 )
  }
  
  if (is.null(xi)) {
    delta <- alpha/sqrt(1+alpha^2)
    xi <- mu - omega * delta * sqrt( 2/pi )
  }
  
  return( list( xi = xi, omega = omega, alpha = alpha ) )
  
}

# ++++++++++++++++++++++++++++++++++++++++++
core_sim <- function( DESIGN, B = 10 ) {  
  
  require(overlapping)
  require(sn)
  # set simulation parameters
  n <- DESIGN$n
  delta <- DESIGN$xi 
  alpha <- DESIGN$alpha 
  omega <- DESIGN$omega
  
  MUSI <- snpar(xi=delta,omega=omega,alpha = alpha)
  
  cat( paste0( "n = ",n,"; mu = ",round(MUSI$mu,1), 
               "; sigma = ",round(MUSI$sigma,2),"; alpha = ",
               round(alpha,2) ),"\n")
  
  # true population overlap
  true_overlap <- integrate( min_dskew_normal, -Inf, Inf, xi = delta, alpha = alpha, 
                             omega = omega )$value
  
  # data simulation
  t( sapply(1:B, function(b){
    
    y1 <- rsn(n,xi=0, omega=1, alpha=0)
    y2 <- rsn(n,xi=delta, omega=omega,alpha=alpha)
    
    mx1 <- mean(y1); sx1 <- sd(y1)
    mx2 <- mean(y2); sx2 <- sd(y2)
    
    ETA1 <- overlap(list(y1,y2))$OV
    ETA2 <- overlap(list(y1,y2), type="2")$OV
    
    norm_pval <- c( shapiro.test(y1)$p.value, shapiro.test(y2)$p.value )
    ks_pval <- ks.test( y1, y2 )$p.value
    vartest_pval <- var.test(y1,y2)$p.value
    wilcox_pval <- wilcox.test(y1,y2)$p.value
    t_pval <- t.test(y1,y2,var.equal = TRUE)$p.value
    welch_pval <- t.test(y1,y2,var.equal = FALSE)$p.value
    
    PERM <- permTest( list(x1=y1,x2=y2), B = NPERM )
    
    output <- c( mx1, sx1, mx2, sx2, ETA1, ETA2, n, delta, alpha, omega, true_overlap, 
                 t_pval, welch_pval, wilcox_pval, vartest_pval,
                 PERM$pval, norm_pval, ks_pval)
    names(output) <- c("mx1", "sx1", "mx2", "sx2", "eta1","eta2", "n", "delta","alpha","omega", 
                       "true_overlap","t_pval","welch_pval",
                       "wilcox_pval","vartest_pval",
                       "zeta_perm_pval","mean_perm_pval","F_perm_pval", 
                       "x1_norm_pval", "x2_norm_pval","ks_test_pval")
    return(output)  
  }) )
}

# ==============================================
# ==============================================
## DESIGN LIST
MUSI <- with( PARlist, expand.grid(mu=mu_vec, sigma=sigma_vec, alpha = alpha_vec) )

DESIGN <- NULL 
for (i in 1:nrow(MUSI)) {
  DESIGN <- rbind( DESIGN, 
    unlist( with( MUSI, sninvpar(mu[i],sigma[i],alpha = alpha[i]) )) )
}

DD <- data.frame(DESIGN)

DESIGN <- NULL
for (n in PARlist$n_vec) {
  DD$n <- n
  DESIGN <- rbind(DESIGN,DD)
}

DESIGN <- apply(DESIGN, 1, as.list)

# ==============================================
# ==============================================
## Simulation
library( parallel )
INIZIO <- Sys.time()
OUT <- mclapply(DESIGN, core_sim, B=B, mc.cores=detectCores())
FINE <- Sys.time()

ALL <- data.frame( do.call(rbind, OUT) )
save(ALL, INIZIO, FINE, B, DESIGN, PARlist, 
     file = paste0(datadir,OUTFILE))

