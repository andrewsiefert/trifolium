
##Parametric bootstrapping for Tweedie GLMM (from cpglmm)
par.boot= function(m, n=1000) {
  y.pred= predict(m, type='response')
  y.rand= replicate(n, rtweedie(length(y.pred), mu=y.pred, power=m$p, phi=m$phi))
  eff.boot= matrix(NA, nrow=n, ncol=length(m$fixef))
  for(i in 1:n) {
    newdat = m@frame
    newdat[,1] = y.rand[,i]
    eff.boot[i,] = update(m, data=newdat)$fixef
  }
  out=data.frame(m$fixef, t(apply(eff.boot, 2, function(i) quantile(i, c(0.025, 0.975)))))
  names(out) = c('coef','lwr','upr')
  return(out)
}

## Plot coefficients and CI's from 2-way interaction
int.2way = function(ci, eff, group) {
  d = ci[grep(paste(":", eff, sep=""), rownames(ci)),]
  d = d + ci[eff,]$coef
  d = rbind(ci[eff,], d)
  d$group = group
  ggplot(d, aes(x=group, y=coef)) + geom_point() + geom_errorbar(aes(ymin=lwr, ymax=upr))
}
  
## Plot coefficients and CI's from 3-way interaction
int.3way = function(ci, eff, g='nrnbr rm', data) {
  effect = ci[eff,]
  nr.eff = ci[paste(g,eff,sep=":"),]
  sp.eff = ci[grep(paste(":", eff, sep=""), rownames(ci)),]
  sp.eff = sp.eff[-(grep(g, rownames(sp.eff))),] 
  sp.nbr.eff = ci[grep(paste(":", g, ":", eff, sep=""), rownames(ci)),]
  
  barb.pres = effect
  sp.pres = sp.eff + effect$coef
  barb.rem = nr.eff + effect$coef
  sp.rem = sp.nbr.eff + sp.pres$coef + nr.eff$coef 
  
  out = rbind(barb.pres, sp.pres, barb.rem, sp.rem)
  out$sp = rep(levels(data$sp), 2)
  out$nr = rep(levels(data$nr), each=length(levels(data$sp)))
  
  dodge = position_dodge(width=0.5)
  ggplot(out, aes(x=sp, y=coef, col=nr)) + geom_point(position=dodge) + 
    geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.2, position=dodge)
}

## Plot fitted vs. observed values
fit.obs= function(mod) {
  y.pred = predict(mod, type='response')
  plot(y.pred, mod@frame[,1]); abline(a=0, b=1); cor(y.pred, mod@frame[,1])^2
}

## Plot observed vs. simulated Tweedie quantiles
qq_tweed <- function(m) {
  
  y_pred <- predict(m, type = 'response')
  r_y <- rtweedie(length(y_pred), mu=y_pred, power=m$p, phi=m$phi)
  y <- m@frame[,1]
  q_obs <- quantile(y, seq(0, 1, 0.01))
  
  d <- list()
  for(i in 1:200) {
    r_y <- rtweedie(length(y_pred), mu=y_pred, power=m$p, phi=m$phi)
    q_tweed <- quantile(r_y, seq(0, 1, 0.01))
    d[[i]] <- tibble('q_obs' = q_obs, 'q_tweed' = q_tweed, rep = i)
  }
  d <- bind_rows(d)
  x_max <- y_max <- max(d$q_obs)
  
  ggplot(d) + geom_line(aes(x = q_obs, y = q_tweed, group = rep), alpha = 0.1, color = "blue") +
    geom_abline(size = 1)
  
}
