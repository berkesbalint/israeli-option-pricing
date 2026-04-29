# =============================================================================
# Izraeli opció árazása kétlépéses Longstaff-Schwartz Monte Carlo módszerrel
# =============================================================================
#
# A függvény izraeli call vagy put opció t = 0-beli árát becsüli a Wang (2024)
# által javasolt kétlépéses LSM eljárással, Black-Scholes részvényár-dinamika
# mellett.
#
# Argumentumok:
#   S0     - kezdő részvényár (S_0)
#   K      - kötési árfolyam
#   r      - kockázatmentes kamatláb (folytonos kamatozás)
#   sigma  - részvényár volatilitása (kappa a dolgozatban)
#   T      - lejárati idő (években)
#   n      - időlépések száma a [0, T] intervallum diszkretizálásához
#   m      - szimulált trajektóriák száma
#   delta  - büntetés, amit az eladó fizet korai lehívás esetén
#   type   - "put" vagy "call" (alapértelmezetten "put")
#   seed   - véletlenszám-generátor magja (opcionális, reprodukálhatósághoz)
#
# Visszatérési érték:
#   Az izraeli opció becsült t = 0-beli ára.
# =============================================================================

price_game_option_two_step_lsmc <- function(S0, K, r, sigma, T, n, m, delta,
                                            type = "put", seed = NULL) {
  
  # --- Bemenetek ellenőrzése ------------------------------------------------
  type <- match.arg(type, choices = c("put", "call"))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  dt <- T / n
  
  # --- Kifizetésfüggvény definiálása az opció típusa szerint ---------------
  # Put: (K - S)^+, Call: (S - K)^+
  payoff <- function(S) {
    if (type == "put") {
      pmax(K - S, 0)
    } else {
      pmax(S - K, 0)
    }
  }
  
  # --- Részvényár-trajektóriák szimulálása ---------------------------------
  # Geometriai Brown-mozgás kockázatsemleges mérték alatt:
  # S_{t+dt} = S_t * exp((r - sigma^2 / 2) * dt + sigma * sqrt(dt) * Z),
  # ahol Z ~ N(0, 1).
  S <- matrix(0, nrow = m, ncol = n + 1)
  S[, 1] <- S0
  
  for (j in 2:(n + 1)) {
    Z <- rnorm(m)
    S[, j] <- S[, j - 1] * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z)
  }
  
  # --- Lejáratkori kifizetések inicializálása ------------------------------
  CF <- payoff(S[, n + 1])
  
  # --- Visszafelé haladó dinamikus programozás -----------------------------
  if (n >= 2) {
    for (t in seq(from = n, to = 2, by = -1)) {
      
      # A következő időlépésből származó kifizetés diszkontált értéke.
      discounted_CF <- exp(-r * dt) * CF
      X <- S[, t]
      
      # Aktuális lehívási érték a vevő (exercise) és az eladó (recall)
      # szempontjából. Az eladó korai lehívásért delta büntetést fizet.
      exercise_payoff <- payoff(X)
      recall_payoff   <- exercise_payoff + delta
      
      # --- 1. lépés: regresszió az OTM trajektóriákon -----------------------
      # Wang (2024) kétlépéses módszerében az OTM pályák közül azokat is
      # bevonjuk a 2. lépés regressziójába, ahol a folytatási érték
      # nagyobb, mint a delta büntetés (azaz az eladónak megérheti várnia).
      N_out_of_money <- which(exercise_payoff == 0)
      N_out_of_sample_relevant <- integer(0)
      
      if (length(N_out_of_money) >= 3) {
        X_otm <- X[N_out_of_money]
        Y_otm <- discounted_CF[N_out_of_money]
        
        # A folytatási érték becslése másodfokú polinomiális regresszióval.
        fit1  <- lm(Y_otm ~ X_otm + I(X_otm^2))
        D_otm <- predict(fit1)
        
        # Csak azokat az OTM pályákat vonjuk be, ahol a folytatás értékesebb,
        # mint a delta (azaz potenciálisan eladói lehívási döntési helyzet).
        N_out_of_sample_relevant <- N_out_of_money[D_otm > delta]
      }
      
      # --- 2. lépés: regresszió a releváns trajektóriákon -------------------
      # ITM pályák + az 1. lépésben kiválasztott OTM pályák.
      N_in_the_money <- which(exercise_payoff > 0)
      N_regression   <- sort(unique(c(N_in_the_money, N_out_of_sample_relevant)))
      
      new_CF <- discounted_CF
      
      if (length(N_regression) >= 3) {
        X_reg <- X[N_regression]
        Y_reg <- discounted_CF[N_regression]
        
        # A folytatási érték (C) becslése a releváns trajektóriákon.
        fit2  <- lm(Y_reg ~ X_reg + I(X_reg^2))
        C_reg <- predict(fit2)
        
        # Vevői lehívás: ha a folytatási érték kisebb, mint az azonnali
        # lehívás kifizetése (a vevő számára megéri lehívni).
        idx_holder_local <- which(C_reg < exercise_payoff[N_regression])
        NExercise_Holder <- N_regression[idx_holder_local]
        
        # Eladói lehívás: ha a folytatási érték nagyobb, mint a recall
        # kifizetés (az eladó számára megéri visszahívni az opciót).
        idx_issuer_local <- which(C_reg > recall_payoff[N_regression])
        NExercise_Issuer <- N_regression[idx_issuer_local]
        
        # A megfelelő trajektóriákon felülírjuk a kifizetéseket.
        new_CF[NExercise_Holder] <- exercise_payoff[NExercise_Holder]
        new_CF[NExercise_Issuer] <- recall_payoff[NExercise_Issuer]
      }
      
      CF <- new_CF
    }
  }
  
  # --- Az opció ára: a t = 1-beli kifizetések diszkontált átlaga ----------
  price <- exp(-r * dt) * mean(CF)
  return(price)
}


# =============================================================================
# Példa-futtatás
# =============================================================================
# price_game_option_two_step_lsmc(S0 = 100, K = 90, r = 0.05, sigma = 0.2,
#                                 T = 1, n = 1000, m = 10000, delta = 10,
#                                 type = "put", seed = 42)
