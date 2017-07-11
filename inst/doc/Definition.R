### R code from vignette source 'Definition.Rnw'

###################################################
### code chunk number 1: changes-and-loss
###################################################
library(data.table)
data(neuroblastoma, package="neuroblastoma")
ids.str <- paste(c(1, 4, 6, 8, 10, 11))
someProfiles <- function(all.profiles){
  data.table(all.profiles)[profile.id %in% ids.str, ]
}
profiles <- someProfiles(neuroblastoma$profiles)
labels <- someProfiles(neuroblastoma$annotations)
problem.list <- split(profiles, profiles[, paste(profile.id, chromosome)])
segs.list <- list()
loss.list <- list()
for(problem.i in seq_along(problem.list)){
  problem.name <- names(problem.list)[[problem.i]]
  pro <- problem.list[[problem.name]]
  meta <- pro[1, .(profile.id, chromosome)]
  max.segments <- min(nrow(pro), 10)
  fit <- tryCatch({
    Segmentor3IsBack::Segmentor(
      pro$logratio, model=2, Kmax=max.segments)
  }, error=function(e){
    print(problem.name)
    NULL
  })
  if(!is.null(fit)){
    rss.vec <- rep(NA, max.segments)
    for(n.segments in 1:max.segments){
      end <- fit@breaks[n.segments, 1:n.segments]
      data.before.change <- end[-n.segments]
      data.after.change <- data.before.change+1
      pos.before.change <- as.integer(
        (pro$position[data.before.change]+pro$position[data.after.change])/2)
      start <- c(1, data.after.change)
      chromStart <- c(pro$position[1], pos.before.change)
      chromEnd <- c(pos.before.change, max(pro$position))
      seg.mean.vec <- fit@parameters[n.segments, 1:n.segments]
      data.mean.vec <- rep(seg.mean.vec, end-start+1)
      residual.vec <- pro$logratio - data.mean.vec
      rss.vec[n.segments] <- sum(residual.vec * residual.vec)
      segs.list[[paste(problem.name, n.segments)]] <- data.table(
        meta,
        n.segments,
        start,
        end,
        chromStart,
        chromEnd,
        mean=seg.mean.vec)
    }
    loss.list[[paste(problem.name, n.segments)]] <- data.table(
      meta,
      n.segments=1:max.segments,
      loss=rss.vec,
      lik=as.numeric(fit@likelihood))
  }
}
loss <- do.call(rbind, loss.list)
segs <- do.call(rbind, segs.list)


###################################################
### code chunk number 2: segs
###################################################
print(segs)


###################################################
### code chunk number 3: loss
###################################################
print(loss)


###################################################
### code chunk number 4: modelSelection
###################################################
selection <- loss[, {
  penaltyLearning::modelSelection(.SD, "loss", "n.segments")
}, by=list(profile.id, chromosome)]
print(selection[profile.id==1 & chromosome==1])


###################################################
### code chunk number 5: modelSelectionplot
###################################################
some.selection <- selection[profile.id==1 & chromosome %in% c(11, 17)]
some.selection[, list(
  pid.chr=paste(profile.id, chromosome),
  min.log.lambda, max.log.lambda, n.segments, loss
)]
library(ggplot2)
gg.selection <- ggplot()+
  theme_grey()+
  facet_grid(. ~ profile.id + chromosome, labeller=label_both)+
  geom_segment(aes(
    min.log.lambda, n.segments,
    xend=max.log.lambda, yend=n.segments),
    data=some.selection)+
  scale_y_continuous(breaks=1:max(some.selection$n.segments))+
  xlab("log(penalty)")
print(gg.selection)


###################################################
### code chunk number 6: changes
###################################################
changes <- segs[1 < start, ]
changes[, list(profile.id, chromosome, n.segments, changepoint=chromStart)]


###################################################
### code chunk number 7: labelError
###################################################
errors <- penaltyLearning::labelError(
  selection, labels, changes,
  change.var="chromStart",
  label.vars=c("min", "max"),
  problem.vars=c("profile.id", "chromosome"))


###################################################
### code chunk number 8: model-errors
###################################################
some.errors <- errors$model.errors[some.selection, on=list(
  profile.id, chromosome, n.segments)]
gg.err <- ggplot()+
  theme_grey()+
  facet_grid(. ~ profile.id + chromosome, labeller=label_both)+
  geom_segment(aes(
    min.log.lambda, errors,
    xend=max.log.lambda, yend=errors),
    data=some.errors)+
  xlab("log(penalty)")
print(gg.err)


###################################################
### code chunk number 9: all-errors
###################################################
all.errors <- data.table(errors$model.errors)
all.errors[, set := ifelse(chromosome=="11", "test", "train")]


###################################################
### code chunk number 10: targetIntervals
###################################################
target.dt <- penaltyLearning::targetIntervals(
  all.errors[set=="train", ],
  c("profile.id", "chromosome"))
target.mat <- target.dt[, cbind(min.log.lambda, max.log.lambda)]
rownames(target.mat) <- target.dt[, paste(profile.id, chromosome)]
print(head(target.mat))


###################################################
### code chunk number 11: featureMatrix
###################################################
feature.dt <- profiles[, list(
  log.data=log(.N),
  log.var=log(median(abs(diff(logratio))))
), by=list(profile.id, chromosome)]
all.feature.mat <- feature.dt[, cbind(log.data, log.var)]
rownames(all.feature.mat) <- feature.dt[, paste(profile.id, chromosome)]
train.feature.mat <- all.feature.mat[rownames(target.mat), ]
print(head(train.feature.mat))


###################################################
### code chunk number 12: IntervalRegression
###################################################
fit <- penaltyLearning::IntervalRegressionUnregularized(
  train.feature.mat, target.mat)
print(fit)


###################################################
### code chunk number 13: pred
###################################################
feature.dt[, pred.log.lambda := predict(fit, all.feature.mat)]


###################################################
### code chunk number 14: ROC
###################################################
test.pred <- feature.dt[chromosome=="11",]
roc <- penaltyLearning::ROChange(
  all.errors, test.pred, c("profile.id", "chromosome"))
pred.thresh <- roc$thresholds[threshold=="predicted"]
gg.roc <- ggplot()+
  geom_path(aes(
    FPR, TPR),
    data=roc$roc)+
  geom_point(aes(
    FPR, TPR),
    data=pred.thresh)
print(gg.roc)
print(pred.thresh)


###################################################
### code chunk number 15: roc-thresh
###################################################
gg.thresh <- ggplot()+
  geom_segment(aes(
    min.thresh, errors,
    xend=max.thresh, yend=errors),
    data=roc$roc)+
  geom_point(aes(
    0, errors),
    data=pred.thresh)+
  xlab("threshold")
print(gg.thresh)


###################################################
### code chunk number 16: vizpred
###################################################
pred.models <- feature.dt[selection, nomatch=0L, on=list(
  profile.id, chromosome,
  pred.log.lambda < max.log.lambda,
  pred.log.lambda > min.log.lambda)]
pred.segs <- segs[pred.models, on=list(profile.id, chromosome, n.segments)]
pred.changes <- pred.segs[1 < start, ]
pred.labels <- errors$label.errors[pred.models, nomatch=0L, on=list(
  profile.id, chromosome, n.segments)]
breakpoint.colors <- c(
  "breakpoint"="#a445ee",
  "normal"="#f6f4bf")
viz.learned <- ggplot()+
  ggtitle("data + labels + learned model segment means and changes")+
  theme_bw()+
  theme(
    legend.position="bottom",
    legend.box="horizontal",
    panel.margin=grid::unit(0, "lines"))+
  facet_grid(profile.id ~ chromosome, scales="free", space="free_x")+
  penaltyLearning::geom_tallrect(aes(
    xmin=min/1e6, xmax=max/1e6, fill=annotation, linetype=status),
    data=pred.labels)+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                                 "false negative"=3,
                                 "false positive"=1))+
  scale_fill_manual("label", values=breakpoint.colors)+
  geom_point(aes(position/1e6, logratio),
             data=profiles,
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))+
  geom_segment(aes(chromStart/1e6, mean, xend=chromEnd/1e6, yend=mean),
               data=pred.segs,
               color="green")+
  geom_vline(aes(xintercept=chromStart/1e6),
             data=pred.changes,
             color="green",
             linetype="dashed")
print(viz.learned)


###################################################
### code chunk number 17: bic-roc
###################################################
bic.dt <- profiles[, list(
  pred.log.lambda=log(log(.N))
), by=list(profile.id, chromosome)]
bic.test <- bic.dt[chromosome==11]
bic.roc <- penaltyLearning::ROChange(
  all.errors, bic.test, c("profile.id", "chromosome"))
print(bic.roc$thresholds[threshold=="predicted"])


###################################################
### code chunk number 18: bic-viz
###################################################
bic.models <- bic.dt[selection, nomatch=0L, on=list(
  profile.id, chromosome,
  pred.log.lambda < max.log.lambda,
  pred.log.lambda > min.log.lambda)]
bic.segs <- segs[bic.models, on=list(profile.id, chromosome, n.segments)]
bic.changes <- bic.segs[1 < start, ]
bic.labels <- errors$label.errors[bic.models, nomatch=0L, on=list(
  profile.id, chromosome, n.segments)]
viz.bic <- ggplot()+
  ggtitle("data + labels + BIC model segment means and changes")+
  theme_bw()+
  theme(
    legend.position="bottom",
    legend.box="horizontal",
    panel.margin=grid::unit(0, "lines"))+
  facet_grid(profile.id ~ chromosome, scales="free", space="free_x")+
  penaltyLearning::geom_tallrect(aes(
    xmin=min/1e6, xmax=max/1e6, fill=annotation, linetype=status),
    data=bic.labels)+
  scale_linetype_manual("error type",
                        values=c(correct=0,
                                 "false negative"=3,
                                 "false positive"=1))+
  scale_fill_manual("label", values=breakpoint.colors)+
  geom_point(aes(position/1e6, logratio),
             data=profiles,
             shape=1)+
  scale_x_continuous(
    "position on chromosome (mega bases)",
    breaks=c(100, 200))+
  geom_segment(aes(chromStart/1e6, mean, xend=chromEnd/1e6, yend=mean),
               data=bic.segs,
               color="green")+
  geom_vline(aes(xintercept=chromStart/1e6),
             data=bic.changes,
             color="green",
             linetype="dashed")
print(viz.bic)


###################################################
### code chunk number 19: scatter
###################################################
scatter.dt <- data.table(
  bic=bic.dt$pred.log.lambda, 
  learned=feature.dt$pred.log.lambda)
gg.scatter <- ggplot()+
  geom_abline(
    slope=1, intercept=0, color="grey")+
  geom_point(aes(
    learned, bic), 
    data=scatter.dt)+
  coord_equal()
print(gg.scatter)


