# =============================================================================
# Izraeli opció árazása binomiális közelítéssel
# =============================================================================
#
# A függvény izraeli call vagy put opció t = 0-beli árát számítja a Cox-Ross-
# Rubinstein-féle binomiális modell segítségével (CRR-közelítés). 
#
# Argumentumok:
#   S0     - kezdő részvényár (S_0)
#   K      - kötési árfolyam
#   r      - kockázatmentes kamatláb (folytonos kamatozás)
#   sigma  - részvényár volatilitása (kappa a dolgozatban)
#   T      - lejárati idő (években)
#   n      - időlépések száma a [0, T] intervallum diszkretizálásához
#   delta  - büntetés, amit az eladó fizet korai lehívás esetén
#   type   - "put" vagy "call" (alapértelmezetten "put")
#
# Visszatérési érték:
#   Az izraeli opció t = 0-beli ára.
# =============================================================================

price_game_option_binomial <- function(S0, K, r, sigma, T, n, delta,
                                       type = "put") {
  
  # --- Bemenetek ellenőrzése ------------------------------------------------
  type <- match.arg(type, choices = c("put", "call"))
  
  # --- A binomiális modell paraméterei (CRR-közelítés) ---------------------
  # A részvényár minden lépésben u-szorosára nő (felfelé) vagy d-szeresére
  # csökken (lefelé), ahol u = exp(sigma * sqrt(dt)) és d = 1 / u.
  # A kockázatsemleges valószínűség: p* = (exp(r*dt) - d) / (u - d).
  dt   <- T / n
  u    <- exp(sigma * sqrt(dt))
  d    <- 1 / u
  disc <- exp(-r * dt)
  p    <- (exp(r * dt) - d) / (u - d)
  
  # --- Részvényár-fa felépítése --------------------------------------------
  # S[i+1, j+1] = a részvény ára az i-edik időpontban, j felfelé lépés után.
  # (R 1-től indexel, ezért a +1 eltolás.)
  S <- matrix(0, n + 1, n + 1)
  for (i in 0:n) {
    for (j in 0:i) {
      S[i + 1, j + 1] <- S0 * u^j * d^(i - j)
    }
  }
  
  # --- Kifizetésfüggvények -------------------------------------------------
  # Y: a vevői (holder) lehívás kifizetése; X = Y + delta: az eladói (issuer)
  # lehívás kifizetése.
  if (type == "put") {
    Y <- pmax(K - S, 0)
  } else {
    Y <- pmax(S - K, 0)
  }
  X <- Y + delta
  
  # --- Az opcióérték fája --------------------------------------------------
  # Lejáratkor (i = n) az opció értéke megegyezik a vevői kifizetéssel,
  # mivel ekkor már nem lehet várni a lehívással.
  V <- matrix(0, n + 1, n + 1)
  V[n + 1, ] <- Y[n + 1, ]
  
  # --- Visszafelé haladó dinamikus programozás -----------------------------
  # Minden korábbi időpontban a Kifer (2000) 3.13. tétel rekurzív formuláját
  # alkalmazzuk:
  #   V_t = min{ X_t, max{ Y_t, disc * E*[V_{t+1}] } },
  # ahol disc * E*[V_{t+1}] = disc * (p * V_up + (1 - p) * V_down) a
  # folytatási érték a kockázatsemleges mérték alatt.
  for (i in (n - 1):0) {
    for (j in 0:i) {
      cont <- disc * (p * V[i + 2, j + 2] + (1 - p) * V[i + 2, j + 1])
      V[i + 1, j + 1] <- min(X[i + 1, j + 1],
                             max(Y[i + 1, j + 1], cont))
    }
  }
  
  # --- Az opció ára t = 0-ban ---------------------------------------------
  return(V[1, 1])
}


# =============================================================================
# Példa-futtatás
# =============================================================================
# price_game_option_binomial(S0 = 100, K = 90, r = 0.05, sigma = 0.2,
#                            T = 1, n = 1000, delta = 10, type = "put")
