plot.divratett.BconsDcons <- function(rates, times, div_only = TRUE, color){
  model <- "BconsDcons"
  i <- which(rates$model == model)
  Bcons <- as.numeric(rates$lambda[i])
  Dcons <- as.numeric(rates$mu[i])
  div <- Bcons-Dcons
  if(div_only == FALSE){
    segments(x0=times[1], x1=times[2], y0=Bcons, y1=Bcons, col= color[model], lwd=1, lty="solid")
    segments(x0=times[1], x1=times[2], y0=Dcons, y1=Dcons, col= color[model], lwd=1, lty="longdash")
  }
  segments(x0=times[1], x1=times[2], y0=div, y1=div, col= color[model], lwd=2, lty="solid")
}

plot.divratett.BexpDcons <- function(rates, times, div_only = TRUE, color){
  model <- "BexpDcons"
  i <- which(rates$model == model)
  Bexp <- function(t) as.numeric(rates$lambda[i]) * exp(as.numeric(rates$alpha[i]) * t)
    Dcons <- as.numeric(rates$mu[i])
    div <- function(t) Bexp(t) - Dcons
    if(div_only == FALSE){
      curve(expr= Bexp, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=1, lty="solid")
      segments(x0=times[1], x1=times[2], y0=Dcons, y1=Dcons, col= color[model], lwd=1, lty="longdash")
    }
    curve(expr= div, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=2, lty="solid")
}

plot.divratett.BlinDcons <- function(rates, times, div_only = TRUE, color){
  model <- "BlinDcons"
  i <- which(rates$model == model)
  Blin <- function(t) as.numeric(rates$lambda[i]) + as.numeric(rates$alpha[i]) * t
  Dcons <- as.numeric(rates$mu[i])
  div <- function(t) Blin(t) - Dcons
  if(div_only == FALSE){
    curve(expr= Blin, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=1, lty="solid")
    segments(x0=times[1], x1=times[2], y0=Dcons, y1=Dcons, col= color[model], lwd=1, lty="longdash")
  }
  curve(expr= div, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=2, lty="solid")
}

plot.divratett.BconsDexp <- function(rates, times, div_only = TRUE, color){
  model <- "BconsDexp"
  i <- which(rates$model == model)
  Bcons <- as.numeric(rates$lambda[i])
  Dexp <- function(t) as.numeric(rates$mu[i]) * exp(as.numeric(rates$beta[i]) * t)
  div <- function(t) Bcons - Dexp(t)
  if(div_only == FALSE){
    segments(x0=times[1], x1=times[2], y0=Bcons, y1=Bcons, col= color[model], lwd=1, lty="solid")
    curve(expr= Dexp, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=1, lty="longdash")
  }
  curve(expr= div, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=2, lty="solid")
}

plot.divratett.BconsDlin <- function(rates, times, div_only = TRUE, color){
  model <- "BconsDlin"
  i <- which(rates$model == model)
  Bcons <- as.numeric(rates$lambda[i])
  Dlin <- function(t) as.numeric(rates$mu[i]) + as.numeric(rates$beta[i]) * t
  div <- function(t) Bcons - Dlin(t)
  if(div_only == FALSE){
    segments(x0=times[1], x1=times[2], y0=Bcons, y1=Bcons, col= color[model], lwd=1, lty="solid")
    curve(expr= Dlin, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=1, lty="longdash")
  }
  curve(expr= div, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=2, lty="solid")
}

plot.divratett.BexpDexp <- function(rates, times, div_only = TRUE, color){
  model <- "BexpDexp"
  i <- which(rates$model == model)
  Bexp <- function(t) as.numeric(rates$lambda[i]) * exp(as.numeric(rates$alpha[i]) * t)
  Dexp <- function(t) as.numeric(rates$mu[i]) * exp(as.numeric(rates$beta[i]) * t)
  div <- function(t) Bexp(t) - Dexp(t)
  if(div_only == FALSE){
    curve(expr= Bexp, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=1, lty="solid")
    curve(expr= Dexp, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=1, lty="longdash")
  }
  curve(expr= div, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=2, lty="solid")
}

plot.divratett.BlinDlin <- function(rates, times, div_only = TRUE, color){
  model <- "BlinDlin"
  i <- which(rates$model == model)
  Blin <- function(t) as.numeric(rates$lambda[i]) + as.numeric(rates$alpha[i]) * t
  Dlin <- function(t) as.numeric(rates$mu[i]) + as.numeric(rates$beta[i]) * t
  div <- function(t) Blin(t) - Dlin(t)
  if(div_only == FALSE){
    curve(expr= Blin, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=1, lty="solid")
    curve(expr= Dlin, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=1, lty="longdash")
  }
  curve(expr= div, from= times[1], to= times[2], add=TRUE, col= color[model], lwd=2, lty="solid")
}
