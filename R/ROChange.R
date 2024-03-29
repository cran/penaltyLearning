ROChange <- structure(function # ROC curve for changepoints
### Compute a Receiver Operating Characteristic curve for a penalty
### function.
(models,
### data.frame describing the number of incorrect labels as a function
### of log(lambda), with columns min.log.lambda, max.log.lambda, fp,
### fn, possible.fp, possible.fn, etc. This can be computed via
### labelError(modelSelection(...), ...)$model.errors -- see examples.
  predictions,
### data.frame with a column named pred.log.lambda, the predicted
### log(penalty) value for each segmentation problem.
  problem.vars=character()
### character: column names used to identify data set / segmentation
### problem.
){
  possible.fp <- possible.fn <- min.log.lambda <- fp <- fn <- thresh <-
    log.lambda <- pred.log.lambda <- errors <- FPR <- tp <- TPR <-
      error.percent <- min.thresh <- max.thresh <- max.log.lambda <-
        next.min <- problems <- n.inconsistent <- min.fp.fn <-
          . <- min.adj <- fp.adj <- fn.adj <- fp.tot.diff <- fn.tot.diff <-
            deriv.min <- deriv.max <-NULL
### The code above is to avoid CRAN NOTEs like
### ROChange: no visible binding for global variable
  if(!(
    is.character(problem.vars) &&
      0 < length(problem.vars) &&
      all(!is.na(problem.vars)) &&
      all(problem.vars %in% names(predictions)) &&
      all(problem.vars %in% names(models))
    )){
    stop("problem.vars should be a character vector of column names (IDs for predictions and models)")
  }
  exp.cols <- c(
    "fp", "possible.fp", "fn", "possible.fn", "errors", "labels",
    problem.vars, "min.log.lambda", "max.log.lambda")
  if(!(
    is.data.frame(models) &&
      all(exp.cols %in% names(models))
    )){
    stop("models should have columns ", paste(exp.cols, collapse=", "))
  }
  for(col.name in exp.cols){
    if(any(is.na(models[[col.name]]))){
      stop(col.name, " should not be NA")
    }
  }
  if(!(
    is.data.frame(predictions) &&
    0 < nrow(predictions) &&
    "pred.log.lambda" %in% names(predictions)
  )){
    stop("predictions should be a data.frame with at least one row and a column named pred.log.lambda")
  }
  if(!all(is.finite(predictions[["pred.log.lambda"]]))){
    stop("all predictions must be finite")
  }
  pred <- data.table(predictions)
  bad.pred <- pred[, list(
    problems=.N
  ), by=problem.vars][1 < problems]
  if(nrow(bad.pred)){
    print(pred[bad.pred, on=problem.vars])
    stop("more than one prediction per problem")
  }
  err <- data.table(models)[pred, on=problem.vars]
  setkeyv(err, c(problem.vars, "min.log.lambda"))
  err.missing <- err[is.na(labels)]
  if(nrow(err.missing)){
    print(err.missing)
    stop("some predictions do not exist in models")
  }
  min.max.not.Inf <- err[, list(
    min=min(min.log.lambda),
    max=max(max.log.lambda)
  ), by=problem.vars][-Inf < min | max < Inf]
  if(nrow(min.max.not.Inf)){
    print(min.max.not.Inf)
    stop("for every problem, the smallest min.log.lambda should be -Inf, and the largest max.log.lambda should be Inf")
  }
  err[, next.min := c(min.log.lambda[-1], Inf), by=problem.vars]
  inconsistent.counts <- err[, list(
    n.inconsistent=sum(next.min != max.log.lambda)
  ), by=problem.vars][0 < n.inconsistent]
  if(nrow(inconsistent.counts)){
    print(err[inconsistent.counts, on=problem.vars])
    stop("max.log.lambda should be equal to the next min.log.lambda")
  }
  for(col.name in c("labels", "possible.fp", "possible.fn")){
    possible.ranges <- err[, {
      x <- .SD[[col.name]]
      list(
        min=min(x),
        max=max(x)
      )}, by=problem.vars]
    possible.inconsistent <- possible.ranges[min != max]
    if(nrow(possible.inconsistent)){
      print(possible.inconsistent)
      stop(
        col.name,
        " should be constant for each problem")
    }
  }
  negative <- err[possible.fp<0 | possible.fn<0 | labels<0]
  if(nrow(negative)){
    print(negative)
    stop("possible.fn/possible.fp/labels should be non-negative")
  }
  possible.name.vec <- c(
    errors="labels",
    fp="possible.fp",
    fn="possible.fn")
  for(err.name in names(possible.name.vec)){
    poss.name <- possible.name.vec[[err.name]]
    poss.num <- err[[poss.name]]
    err.num <- err[[err.name]]
    out.of.range <- err[poss.num < err.num | err.num < 0]
    if(nrow(out.of.range)){
      print(out.of.range)
      stop(
        err.name,
        " should be in [0,",
        poss.name,
        "]")
    }
  }
  first.dt <- err[min.log.lambda==-Inf]
  last.dt <- err[max.log.lambda==Inf]
  total.dt <- first.dt[, list(
    labels=sum(labels),
    possible.fp=sum(possible.fp),
    possible.fn=sum(possible.fn))]
  if(total.dt$possible.fp==0){
    stop("no negative labels")
  }
  if(total.dt$possible.fn==0){
    stop("no positive labels")
  }
  thresh.dt <- err[order(min.log.lambda), {
    fp.prb.diff <- diff(fp)
    fn.prb.diff <- diff(fn)
    any.change <- fp.prb.diff != 0 | fn.prb.diff != 0
    data.table(
      log.lambda=max.log.lambda[c(any.change, FALSE)],
      fp.prb.diff=as.numeric(fp.prb.diff[any.change]),
      fn.prb.diff=as.numeric(fn.prb.diff[any.change]))
  }, by=problem.vars]
  fp.fn.problems <- thresh.dt[pred, on=problem.vars, nomatch=0L]
  fp.fn.problems[, thresh := log.lambda - pred.log.lambda]
  uniq.thresh <- fp.fn.problems[order(thresh), list(
    fp.tot.diff=sum(fp.prb.diff),
    fn.tot.diff=sum(fn.prb.diff)
  ), by=thresh]
  fp.fn.totals <- uniq.thresh[, data.table(
    min.thresh=c(-Inf, thresh),
    max.thresh=c(thresh, Inf),
    fp = cumsum(c(sum(first.dt$fp), fp.tot.diff)),
    fn = rev(cumsum(rev(c(-fn.tot.diff, sum(last.dt$fn)))))
  )]
  ## Compute aum = area under min(fp,fn).
  fp.fn.totals[, min.fp.fn := pmin(fp, fn)]
  aum <- fp.fn.totals[, sum(ifelse(
    min.fp.fn==0, 0, min.fp.fn*(max.thresh-min.thresh)))]
  ## Compute directional derivatives coming from lo and hi values. The
  ## main idea is that we look at what happens if the predicted value
  ## for a particular problem is decreased, which results in a bigger
  ## threshold. If that threshold/diff is relevant then it will result
  ## in a change in the min after the threshold, relative to the
  ## actual min after the threshold.
  pos.or.neg.vec <- c(min=-1, max=1)
  for(min.or.max in names(pos.or.neg.vec)){
    pos.or.neg <- pos.or.neg.vec[[min.or.max]]
    ## To compute the directional derivatives of aum we need to join
    ## fp.fn.totals(min.thresh,max.thresh,fp,fn) to the diffs with
    ## individual thresholds/problems, fp.fn.problems(
    ## problem,log.lambda,fp.prb.diff,fn.prb.diff,pred.log.lambda,thresh)
    fp.fn.join <- fp.fn.totals[
      fp.fn.problems,
      .(fp, fn, min.fp.fn, fp.prb.diff, fn.prb.diff),
      on=structure("thresh", names=paste0(min.or.max, ".thresh"))]
    ## changes fp.prb.diff,fn.prb.diff occur at thresh.
    ## fp,fn are the totals in (min.thresh,max.thresh).
    ## by joining min.thresh or max.thresh = thresh,
    ## we get the total fp/fn after or before thresh.
    for(fX in c("fp", "fn")){
      ## here we compute what happens to the fp/fn adjacent to the
      ## totals, on the other side of thresh, when we apply
      ## fp.prb.diff,fn.prb.diff for each thresh/problem.  e.g. when
      ## min.or.max=="min" the fp/fn columns give the total fp/fn
      ## after the thresh where there is a change in the problem fp/fn
      ## curve, and so here we compute the fp/fn before the thresh.
      set(
        fp.fn.join,
        j=paste0(fX, ".adj"),
        value=fp.fn.join[[fX]]+
          pos.or.neg*fp.fn.join[[paste0(fX, ".prb.diff")]]
      )
    }
    ## here we compute Min(FP,FN) adjacent to the totals, on the other
    ## side of thresh, when we apply fp.prb.diff,fn.prb.diff for each
    ## thresh/problem.
    fp.fn.join[, min.adj := pmin(fp.adj, fn.adj)]
    ## finally we check to see how much Min(FP,FN) would change via
    ## fp.prb.diff,fn.prb.diff at each thresh/problem. this amounts to
    ## computing min.after-min.before the thresh.
    set(
      fp.fn.problems,
      j=paste0("deriv.", min.or.max),
      value=fp.fn.join[, pos.or.neg*(min.adj - min.fp.fn)]
    )
  }
  ## Compute TPR/FPR rates for ROC-AUC analysis.
  interval.dt <- fp.fn.totals[, data.table(
    total.dt,
    min.thresh,
    max.thresh,
    fp, fn, min.fp.fn)]
  interval.dt[, errors := fp+fn]
  interval.dt[, FPR := fp/possible.fp]
  interval.dt[, tp := possible.fn - fn]
  interval.dt[, TPR := tp/possible.fn]
  interval.dt[, error.percent := 100*errors/labels]
  dist00.vec <- interval.dt[c(1, .N), sqrt(FPR^2+TPR^2)]
  indices <- if(dist00.vec[1]<dist00.vec[2]){
    1:nrow(interval.dt)
  }else{
    nrow(interval.dt):1
  }
  sorted.dt <- interval.dt[indices]
  roc.polygon <- sorted.dt[, {
    has11 <- FPR[.N]==1 & TPR[.N]==1
    has00 <- FPR[1]==0 & TPR[1]==0
    list(
      FPR=c(if(!has00)0, FPR, if(!has11)1, 1),
      TPR=c(if(!has00)0, TPR, if(!has11)1, 0)
      )
  }]
  # if this is a sequence from q=1 to Q subscripts are:
  left <- roc.polygon[-.N] # _q
  right <- roc.polygon[-1] # _{q+1}
  ##value<< named list of results:
  list(
    roc=interval.dt, ##<< a data.table with one row for each point on
                     ##the ROC curve
    thresholds=rbind(##<< two rows of roc which correspond to the
                     ##predicted and minimal error thresholds
      data.table(
        threshold="predicted",
        interval.dt[min.thresh < 0 & 0 <= max.thresh, ]),
      data.table(threshold="min.error", interval.dt[which.min(errors), ])),
    auc.polygon=roc.polygon, ##<<a data.table with one row for
                             ##each vertex of the polygon used to
                             ##compute AUC
    auc=##<<numeric Area Under the ROC curve
      sum((right$FPR-left$FPR)*(right$TPR+left$TPR)/2),
    aum=aum, ##<< numeric Area Under Min(FP,FN)
    aum.grad=##<< data.table with one row for each prediction, and
      ##columns hi/lo bound for the aum
      ##generalized gradient.
      fp.fn.problems[pred, .(
        lo=sum(deriv.min, na.rm=TRUE),
        hi=sum(deriv.max, na.rm=TRUE)
      ), by=.EACHI, on=problem.vars]
  )
  ##end<<
}, ex=function(){

  library(penaltyLearning)
  library(data.table)

  data(neuroblastomaProcessed, envir=environment())
  ## Get incorrect labels data for one profile.
  pid <- 11
  pro.errors <- neuroblastomaProcessed$errors[
    profile.id==pid][order(chromosome, min.log.lambda)]
  dcast(pro.errors, n.segments ~ chromosome, value.var="errors")
  ## Get the feature that corresponds to the BIC penalty = log(n),
  ## meaning log(penalty) = log(log(n)).
  chr.vec <- paste(c(1:4, 11, 17))
  pid.names <- paste0(pid, ".", chr.vec)
  BIC.feature <- neuroblastomaProcessed$feature.mat[pid.names, "log2.n"]
  pred <- data.table(pred.log.lambda=BIC.feature, chromosome=chr.vec)
  ## edit one prediction so that it ends up having the same threshold
  ## as another one, to illustrate an aum sub-differential with
  ## un-equal lo/hi bounds.
  err.changes <- pro.errors[, {
    .SD[c(NA, diff(errors) != 0), .(min.log.lambda)]
  }, by=chromosome]
  (ch.vec <- err.changes[, structure(min.log.lambda, names=chromosome)])
  other <- "11"
  (diff.other <- ch.vec[[other]]-pred[other, pred.log.lambda, on=.(chromosome)])
  pred["1", pred.log.lambda := ch.vec[["1"]]-diff.other, on=.(chromosome)]
  pred["4", pred.log.lambda := 2, on=.(chromosome)]
  ch.vec[["1"]]-pred["1", pred.log.lambda, on=.(chromosome)]
  result <- ROChange(pro.errors, pred, "chromosome")
  library(ggplot2)
  ## Plot the ROC curves.
  ggplot()+
    geom_path(aes(FPR, TPR), data=result$roc)+
    geom_point(aes(FPR, TPR, color=threshold), data=result$thresholds, shape=1)

  ## Plot the number of incorrect labels as a function of threshold.
  ggplot()+
    geom_segment(aes(
      min.thresh, errors,
      xend=max.thresh, yend=errors),
      data=result$roc)+
    geom_point(aes((min.thresh+max.thresh)/2, errors, color=threshold),
               data=result$thresholds,
               shape=1)+
    xlab("log(penalty) constant added to BIC penalty")

  ## Plot area under Min(FP,FN).
  err.colors <- c(
    "fp"="red",
    "fn"="deepskyblue",
    "min.fp.fn"="black")
  err.sizes <- c(
    "fp"=3,
    "fn"=2,
    "min.fp.fn"=1)
  roc.tall <- melt(result$roc, measure.vars=names(err.colors))
  area.rects <- data.table(
    chromosome="total",
    result$roc[0<min.fp.fn])
  (gg.total <- ggplot()+
     geom_vline(
       xintercept=0,
       color="grey")+
     geom_rect(aes(
       xmin=min.thresh, xmax=max.thresh,
       ymin=0, ymax=min.fp.fn),
       data=area.rects,
       alpha=0.5)+
     geom_text(aes(
       min.thresh, min.fp.fn/2,
       label=sprintf(
         "Area Under Min(FP,FN)=%.3f ",
         result$aum)),
       data=area.rects[1],
       hjust=1,
       color="grey50")+
     geom_segment(aes(
       min.thresh, value,
       xend=max.thresh, yend=value,
       color=variable, size=variable),
       data=data.table(chromosome="total", roc.tall))+
     scale_size_manual(values=err.sizes)+
     scale_color_manual(values=err.colors)+
     theme_bw()+
     theme(panel.grid.minor=element_blank())+
     scale_x_continuous(
       "Prediction threshold")+
     scale_y_continuous(
       "Incorrectly predicted labels",
       breaks=0:10))

  ## Add individual error curves.
  tall.errors <- melt(
    pro.errors[pred, on=.(chromosome)],
    measure.vars=c("fp", "fn"))
  gg.total+
    geom_segment(aes(
      min.log.lambda-pred.log.lambda, value,
      xend=max.log.lambda-pred.log.lambda, yend=value,
      size=variable, color=variable),
      data=tall.errors)+
    facet_grid(chromosome ~ ., scales="free", space="free")+
    theme(panel.spacing=grid::unit(0, "lines"))+
    geom_blank(aes(
      0, errors),
      data=data.table(errors=c(1.5, -0.5)))

  print(result$aum.grad)
  if(interactive()){#this can be too long for CRAN.
    ## Plot how Area Under Min(FP,FN) changes with each predicted value.
    aum.dt <- pred[, {
      data.table(log.pen=seq(0, 4, by=0.5))[, {
        chr <- paste(chromosome)
        new.pred.dt <- data.table(pred)
        new.pred.dt[chr, pred.log.lambda := log.pen, on=.(chromosome)]
        with(
          ROChange(pro.errors, new.pred.dt, "chromosome"),
          data.table(aum))
      }, by=log.pen]
    }, by=chromosome]
    bounds.dt <- melt(
      result$aum.grad,
      measure.vars=c("lo", "hi"),
      variable.name="bound",
      value.name="slope")[pred, on=.(chromosome)]
    bounds.dt[, intercept := result$aum-slope*pred.log.lambda]
    ggplot()+
      geom_abline(aes(
        slope=slope, intercept=intercept),
        size=1,
        data=bounds.dt)+
      geom_text(aes(
        2, 2, label=sprintf("directional derivatives = [%d, %d]", lo, hi)),
        data=result$aum.grad)+
      scale_color_manual(
        values=c(
          predicted="red",
          new="black"))+
      geom_point(aes(
        log.pen, aum, color=type),
        data=data.table(type="new", aum.dt))+
      geom_point(aes(
        pred.log.lambda, result$aum, color=type),
        shape=1,
        data=data.table(type="predicted", pred))+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      facet_wrap("chromosome", labeller=label_both)+
      coord_equal()+
      xlab("New log(penalty) value for chromosome")+
      ylab("Area Under Min(FP,FN)
using new log(penalty) for this chromosome
and predicted log(penalty) for others")
  }

})

