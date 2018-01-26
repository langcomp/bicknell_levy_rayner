################################################################################
### etAnalyze.R
###
### created by Roger Levy
### additionally developed by Klinton Bicknell and Emily Morgan
################################################################################
library(parallel)
library('methods')

sentence.end <- 9999                    # upper bound on line length

setClass("fixations",representation(x="numeric",
                                    y="numeric",
                                    region="numeric",
                                    pos="numeric", # position within region
                                    startTime="numeric",
                                    endTime="numeric"))
validFixations <- function(object) {
  len <- length(object@x)
  if(length(object@y) != len ||
     length(object@startTime) != len ||
     length(object@endTime) != len)
    return("mismatch in length of slots!")
  else  return(TRUE)
}
setValidity("fixations",validFixations)
setClass("trial",representation(condition="numeric",
                                order="numeric",
                                subj="numeric",
                                item="numeric",
                                totalTime="numeric"), contains="fixations")
validTrial <- function(object) {
  r <- validFixations(object)
  if(r != TRUE)
    return(r)
  if(length(object@subj) != 1 ||
     length(object@item) != 1 ||
     length(object@condition) != 1 ||
     length(object@order) != 1 ||
     length(object@totalTime) != 1)
    return("trial-global variable too long!!")
  else  return(TRUE)
}
setValidity("trial", validTrial)

setClass("ETexpt", representation(trials="list", delim="array",
                                  line.length="numeric"))

## This function is now deprecated. Instead, use readTrials below, which is
## vectorized and has new functionality. (as of 28 Feb 2009, only works for
## single-line experiments.)
readTrialsFromFile <- function(filename, subj.number, delim, line.length,
                               exclude.unknown.fixations=TRUE,
                               exclude.conditions=c(), include.items=c(),
                               min.fix=80) {
  n <- length(readLines(filename))
  result <- list()
  for(i in 1:n) {
    line <- scan(filename,skip=i-1,nlines=1,quiet=T)
    order <- line[1]
    cond <- line[2]
    if(cond %in% exclude.conditions)
      next
    item <- line[3]
    if (length(include.items)>0 & !(item %in% include.items))
    	next
    tot.time <- line[4]
    nfix <- line[8]
    fixs <- line[9:length(line)]
    if(length(fixs) %% 4 != 0)
      return(paste("error -- fixations vector wrong length on line!", i))
    ##cat(nfix,"\n")
    ##cat(length(fixs)/ 4,"\n")
    x <- fixs[(1:nfix)*4 - 3]
    y <- fixs[(1:nfix)*4 - 2]
    startTime <- fixs[(1:nfix)*4 - 1]
    endTime <- fixs[(1:nfix)*4]
    fixTime <- endTime-startTime
    ok.fix <- fixTime >= min.fix
    x <- x[ok.fix]
    y <- y[ok.fix]
    startTime <- startTime[ok.fix]
    endTime <- endTime[ok.fix]
    if(exclude.unknown.fixations) {
      idx <- x>=0 & y >= 0
      x <- x[idx]
      y <- y[idx]
      startTime <- startTime[idx]
      endTime <- endTime[idx]
    }
    this.delim <- delim[item,cond,1:dim(delim)[3]]
    this.region <- sapply(x + y*line.length,function(z) match(TRUE, z < c(this.delim,1000)) - 1)   ## we'd need to change this in order to deal with multi-line experiments
    this.pos <- x - delim[this.region]
    
    if(length(this.region)==0)
      next
    result <- append(result,new("trial",condition=cond,order=order,subj=subj.number, item=item,totalTime=tot.time,x=x,y=y,region=this.region,startTime=startTime,endTime=endTime,pos=this.pos))
  }
  return(result)
}

## this reads a single da1 file. (it may not work for multi-line experiments)
readTrials <- function(filename, subj.number, delim, line.length,
                       exclude.unknown.fixations=TRUE,
                       ## unknown fixations are cases where the x or y values
                       ## are < 0. I'm [KB] guessing these correspond to cases
                       ## in which the tracker places the fixation off-screen?
                       exclude.conditions=c(), include.items=c(),
                       min.fix=80,      # in ms
                       ignore.num.fixes=FALSE # needed because timdrop.pl can
                                               # break the nfixes column
                       ) {
  file.list <- lapply(strsplit(readLines(filename), '\\s+'), as.integer)
  headings <- sapply(file.list, '[', c(1:4,8))
  orders <- headings[1,]
  conds <- headings[2,]
  items <- headings[3,]
  tot.times <- headings[4,]
  nfixes <- headings[5,]
  fixes <- lapply(file.list,
                  function(x) {
                    len <- length(x)
                    if (len >= 9) {
                      return(x[9:len])
                    } else {
                      return(numeric())
                    }
                  })
    
  if (!ignore.num.fixes) {
    stopifnot(isTRUE(all.equal(nfixes, sapply(fixes, length)/4)))
  }
  fixations <- lapply(fixes,
                      function(x) {
                        stopifnot(length(x) %% 4 == 0)
                        x <- data.frame(matrix(x,ncol=4,byrow=T))
                        names(x) <- c('x','y','startTime','endTime')
                        x$duration <- x$endTime-x$startTime
                        return(x)
                      })

  if (length(exclude.conditions) > 0) {
    to.exclude <- conds %in% exclude.conditions
    fixations <- fixations[!to.exclude]
    orders <- orders[!to.exclude]
    conds <- conds[!to.exclude]
    items <- items[!to.exclude]
    tot.times <- tot.times[!to.exclude]
  }

  if (length(include.items) > 0) {
    to.include <- items %in% include.items
    fixations <- fixations[to.include]
    orders <- orders[to.include]
    conds <- conds[to.include]
    items <- items[to.include]
    tot.times <- tot.times[to.include]
  }

  ## Do all relevant exclusions on individual fixations here. This way, we can
  ## exclude trials that have no valid fixations left.
  if (exclude.unknown.fixations) {
    fixations <- lapply(fixations,
                        function(f) {
                          return(subset(f, x >= 0 & y >= 0))
                        })
  }
  if (min.fix > 0) {
    fixations <- lapply(fixations,
                        function(f) {
                          return(subset(f, duration >= min.fix))
                        })
  }

  ## now get rid of trials with no fixations
  nonempty.trials <- sapply(fixations, function(f) nrow(f)>0)
  fixations <- fixations[nonempty.trials]
  orders <- orders[nonempty.trials]
  conds <- conds[nonempty.trials]
  items <- items[nonempty.trials]
  tot.times <- tot.times[nonempty.trials]

  trials <- mapply(function(fixations,order,cond,item,tot.time) {
    this.delim <- delim[item, cond, ]
    fixations$region <-
      mapply(
        function(x,y) {
          match(T, x+y*line.length < c(delim[item,cond,],sentence.end))-1
        }, fixations$x, fixations$y)
    fixations$pos <- fixations$x - this.delim[fixations$region]

    return(new("trial", condition=cond, order=order, subj=subj.number,
               item=item, totalTime=tot.time, x=fixations$x, y=fixations$y,
               region=fixations$region, startTime=fixations$startTime,
               endTime=fixations$endTime, pos=fixations$pos))
  }, fixations, orders, conds, items, tot.times)
  
  return(trials)
}

readDelim <- function(filename) {
  ## assumes item numbers and cond numbers start at 1
  file.list <- readLines(filename)
  delim.list <- lapply(strsplit(file.list, '\\s+'), as.integer)
  ## get rid of initial NA caused by lines beginning with space
  delim.list <- lapply(delim.list,
                     function(x) {
                       if (is.na(x[1])) {
                         x <- x[2:length(x)]
                       }
                       return(x)
                     })
  headings <- sapply(delim.list, '[', 1:3)
  items <- headings[1,]
  conds <- headings[2,]
  num.regions <- headings[3,]
  max.region <- max(num.regions)
  max.cond <- max(conds)
  max.item <- max(items)
  stopifnot(isTRUE(all.equal(num.regions,
                             unlist(lapply(delim.list, length))-3)),
            max.cond*max.item == length(delim.list))
  delim <- t(sapply(delim.list,
                    function(x) {x <- c(x,rep(sentence.end,
                                              max.region-(length(x)-3)))
                                 return(x)}))
  delim <- delim[order(delim[, 1], delim[, 2]), 4:ncol(delim)]
  dim(delim) <- c(max.cond, max.item, max.region)
  delim <- aperm(delim, c(2,1,3)) # now [item, cond, region]
  return(delim)
}

getETexpt <- function(da1.files, subj.numbers, delim.file, line.length=160,
                      exclude.conditions=c(), include.items=c(), min.fix=80,
                      ignore.num.fixes=F,
                      mc.cores=getOption("mc.cores", 2L) # use 1 for no forking
                      ) {
  #cat("Getting ET experiment from files", da1.files, "...")
  delim <- readDelim(delim.file)
  stopifnot(length(da1.files) == length(subj.numbers))

  trials <- mcmapply(
    function(da1, subj) {
      readTrials(da1, subj, delim, line.length,
                 exclude.conditions=exclude.conditions,
                 include.items=include.items,
                 min.fix=min.fix,
                 ignore.num.fixes=ignore.num.fixes)
    }, da1.files, subj.numbers, SIMPLIFY=F, USE.NAMES=F, mc.cores=mc.cores)
  
  trials <- do.call(c, trials)

  return(new("ETexpt",trials=trials,delim=delim,line.length=line.length))
}

standardize.df <- function(df) {
  ### puts data frames in a common format to allow for easy merging of different
  ### ET measures.
  df <- df[do.call(order, df), ] # sort by all columns
  row.names(df) <- NULL
  return(df)
}

## total.time <- function(expt, region) {
##  
## }

all.saccade.lengths <- function(expt) {
  subj <- lapply(expt@trials, function(x) {x@subj})
  item <- lapply(expt@trials, function(x) {x@item})
  cond <- lapply(expt@trials, function(x) {x@condition})
  ord <- lapply(expt@trials, function(x) {x@order})
  saclen <- lapply(expt@trials, function(t) {
    saclens <- c(t@x[2:length(t@x)],NA) - t@x
    return(saclens[1:length(saclens)-1])
  })
  result <- mapply(
    function(subj, item, order, cond, saclen) {
      data.frame(subj,item,order,cond,saclen)
    }, subj, item, cond, ord, saclen, SIMPLIFY=F)
  result <- do.call(rbind, result)
  return(standardize.df(result))
}

landing.positions <- function(expt, region) {
  subj <- sapply(expt@trials, function(x) {x@subj})
  item <- sapply(expt@trials, function(x) {x@item})
  cond <- sapply(expt@trials, function(x) {x@condition})
  ord <- sapply(expt@trials, function(x) {x@order})
  lpos <- sapply(expt@trials, function(t) {
    max.region <- 0
    for(j in 1:length(t@x)) {
      if(max.region < t@region[j] & t@region[j] == region) {
        return(t@pos[j])
      }
      max.region <- max(max.region,t@region[j])
    }
    return(-1)
  })
  result <- data.frame(subj,item,order=ord,cond,lpos)
  return(standardize.df(result))
}

## gives info about the first fixation made on a region and the
## saccade that was launched from there (regardless of where it went)
second.saccade.info <- function(expt, region) {
  subj <- sapply(expt@trials, function(x) {x@subj})
  item <- sapply(expt@trials, function(x) {x@item})
  cond <- sapply(expt@trials, function(x) {x@condition})
  ord <- sapply(expt@trials, function(x) {x@order})
  results <- lapply(expt@trials, function(t) {
    max.region <- 0
    lpos <- NA
    x <- NA
    x.next <- NA
    region.next <- NA
    lpos.next <- NA
    for(j in 1:length(t@x)) {
      if(max.region < t@region[j] & t@region[j] == region) {
        lpos <- t@pos[j]
        x <- t@x[j]
        if (j != length(t@x)) {
            x.next <- t@x[j+1]
            region.next <- t@region[j+1]
            lpos.next <- t@pos[j+1]
        }
        break
      }
      max.region <- max(max.region,t@region[j])
    }
    return(data.frame(lpos=lpos, x=x, x.next=x.next,
                      region.next=region.next,
                      lpos.next=lpos.next))
  })
  results <- do.call(rbind, results)
  result <- data.frame(subj, item, order=ord, cond)
  result <- cbind(result, results)
  return(standardize.df(result))
}

launch <- function(expt, region) {
  subj <- sapply(expt@trials, function(x) {x@subj})
  item <- sapply(expt@trials, function(x) {x@item})
  cond <- sapply(expt@trials, function(x) {x@condition})
  ord <- sapply(expt@trials, function(x) {x@order})
  launch <- sapply(expt@trials, function(t) {
    max.region <- 0
    for(j in 1:length(t@x)) {
      if(max.region < t@region[j] & t@region[j] == region) {
        if (j>1) {
          saclen <- t@x[j]-t@x[j-1]
          lpos <- t@pos[j]
          return(lpos-saclen)
        } else {
          return(NA)
        }
      }
      max.region <- max(max.region,t@region[j])
    }
    return(NA)
  })
  result <- data.frame(subj,item,order=ord,cond,launch)
  return(standardize.df(result))
}

## returns a data frame of first fixations
first.fixation.durations <- function(expt, region) {
  subj <- c()
  item <- c()
  cond <- c()
  ord <- c()
  rt <- c()
  for(i in 1:length(expt@trials)) {
    if(region > dim(expt@delim)[3]+1)
      cat("Error -- region number ", region, " too high.\n")
    t <- expt@trials[[i]]
    max.region <- 0
    this.rt <- 0
    for(j in 1:length(t@x)) {
      if(max.region < t@region[j] & t@region[j] == region) {
        this.rt <- t@endTime[j] - t@startTime[j]
        break
      }
      max.region <- max(max.region,t@region[j])
    }
    subj <- c(subj, t@subj)
    item <- c(item, t@item)
    cond <- c(cond, t@condition)
    ord <- c(ord, t@order)
    rt <-  c(rt,this.rt)
  }
  result <- data.frame(subj=subj,item=item,order=ord,cond=cond,rt=rt)
  return(standardize.df(result))
}


## some clear cases of mismatch thus far...(14 Dec 08, RPL)
## region "3" (2 in eyedry)
## 1,11; 1,18
first.pass.durations <- function(expt,region) {
  subj <- c()
  item <- c()
  cond <- c()
  ord <- c()
  rt <- c()
  for(i in 1:length(expt@trials)) {
    if(region > dim(expt@delim)[3])
      warning("Error -- region number ", region, " too high.\n")
    t <- expt@trials[[i]]
    subj <- c(subj, t@subj)
    item <- c(item, t@item)
    cond <- c(cond, t@condition)
    ord <- c(ord, t@order)
    leftEdge <- expt@delim[t@item, t@condition, region]
    rightEdge <- ifelse(region==(dim(expt@delim)[3]), 999,
                        expt@delim[t@item,t@condition, region+1])
    passedLeftEdge <- t@x >= leftEdge
    passedRightEdge <- t@x >= rightEdge
    this.rt <- 0
    for(j in 1:length(t@x)) {
      if(passedRightEdge[j])
        break
      if((! passedLeftEdge[j]) & this.rt > 0)
        break
      if(passedLeftEdge[j])
        this.rt <- this.rt + t@endTime[j] - t@startTime[j]
    }
    rt <- c(rt,this.rt)
  }
  result <- data.frame(subj=subj,item=item,order=ord,cond=cond,rt=rt)
  return(standardize.df(result))
}


### Still needs to be stress-test compared against eyedry output
gopast.durations <- function(expt,region) {
  subj <- c()
  item <- c()
  cond <- c()
  ord <- c()
  rt <- c()
  for(i in 1:length(expt@trials)) {
    if(region > dim(expt@delim)[3]+1)
      warning("Error -- region number ", region, " too high.\n")
    t <- expt@trials[[i]]
    subj <- c(subj, t@subj)
    item <- c(item, t@item)
    cond <- c(cond, t@condition)
    ord <- c(ord, t@order)
    leftEdge <-  expt@delim[t@item, t@condition, region]
    rightEdge <- ifelse(region==(dim(expt@delim)[3]), 999,
                        expt@delim[t@item,t@condition, region+1])
    passedLeftEdge <- t@x >= leftEdge
    passedRightEdge <- t@x >= rightEdge
    previouslyPassedLeftEdge <- FALSE
    this.rt <- 0
    for(j in 1:length(t@x)) {
      if(passedLeftEdge[j]) {
        previouslyPassedLeftEdge <- TRUE
      }
      if(passedRightEdge[j])
        break
      if(previouslyPassedLeftEdge)
         this.rt <- this.rt + t@endTime[j] - t@startTime[j]
    }
    rt <- c(rt,this.rt)
  }
  result <- data.frame(subj=subj,item=item,order=ord,cond=cond,rt=rt)
  return(standardize.df(result))
}

## pretty close to matching eyedry output -- 5 out of 819 cases mismatch (8 Sep
## 09 -- RPL)
first.pass.regressions <- function(expt) {
  regions <- 1:dim(expt@delim)[3]
  subj <- c()
  item <- c()
  cond <- c()
  reg <- c()
  ord <- c()
  fpr <- c()
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    previous.region <- 0
    this.trial.regressions <- rep(FALSE,length(regions))
    first.pass.completed <- rep(FALSE,length(regions))
    current.first.pass.region <- 0
    for(j in 1:length(t@x)) {
      if(t@region[j] < previous.region && ! first.pass.completed[previous.region])
        this.trial.regressions[previous.region] <- TRUE
      if(t@region[j] != previous.region)
        for(z in 1:previous.region)
          first.pass.completed[z] <- TRUE
      previous.region <- t@region[j]
    }
    for(j in 1:length(this.trial.regressions)) {
      subj <- c(subj,t@subj)
      item <- c(item,t@item)
      cond <- c(cond,t@condition)
      ord <- c(ord,t@order)
      reg <- c(reg,j)
      fpr <- c(fpr,this.trial.regressions[j])
    }
  }
  result <- data.frame(subj=subj,item=item,order=ord,cond=cond,region=reg,fpr=fpr)
  return(standardize.df(result))
}

second.pass.rt <- function(expt, skip.counts.as.first.pass=TRUE) {
  regions <- 1:dim(expt@delim)[3]
  subj <- numeric(0)
  item <- numeric(0)
  cond <- numeric(0)
  reg <- numeric(0)
  ord <- numeric(0)
  second.pass <- c()
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    previous.region <- 0
    this.trial.rt <- rep(0,length(regions))
    first.pass.completed <- rep(FALSE,length(regions))
    visited <- rep(FALSE,length(regions))
    for(j in 1:length(t@x)) {
      if(t@region[j] < 1) {
        print("error: trial with region 0: ")
        print(t)
      }
      visited[t@region[j]] <- TRUE
# old version:
#      if(t@region[j] != previous.region)
#        first.pass.completed[previous.region] <- TRUE
#      if(skip.counts.as.first.pass)
#        first.pass.completed[1:(t@region[j]-1)] <- TRUE
#      if(first.pass.completed[t@region[j]])
#        this.trial.rt[t@region[j]] <- this.trial.rt[t@region[j]] + t@endTime[j] - t@startTime[j]
      if(t@region[j] > previous.region) {
        #first.pass.completed[previous.region] <- TRUE
        if(skip.counts.as.first.pass)
        	first.pass.completed[1:(t@region[j]-1)] <- TRUE
        else { #first.pass.completed <- TRUE for all the regions to the left that you've visited
        	for (k in 1:t@region[j]-1) {
        		if (visited[k])
        			first.pass.completed[k] <- TRUE
        	}
       	}
      }
      if(first.pass.completed[t@region[j]])
        this.trial.rt[t@region[j]] <- this.trial.rt[t@region[j]] + t@endTime[j] - t@startTime[j]

      previous.region <- t@region[j]
    }
    for(j in 1:length(this.trial.rt)) {
      subj <- c(subj,t@subj)
      item <- c(item,t@item)
      cond <- c(cond,t@condition)
      ord <- c(ord,t@order)
      reg <- c(reg,j)
      second.pass <- c(second.pass,this.trial.rt[j])
    }
  }
  result <- data.frame(subj=subj,item=item,order=ord,cond=cond,region=reg,second.pass=second.pass)
  return(standardize.df(result))
}

qa.accuracy <- function(expt) {
  do.call('rbind', lapply(expt@trials, function(x) list()))
}

## computes relative frequencies of first-pass regressions from region
## \code{region} to all previous regions type ranges over count, prop,
## posterior.ci.lower, posterior.ci.upper (latter 2 with uniform beta-binomial
## model)
regression.target.breakdown <- function(expt, conds, originating.region,
                                        type="count") {
  result <- numeric(dim(expt@delim)[3])
  total.regressions <- 0
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    if(! t@condition %in% conds)
      next;
    highest.region <- 0
    previous.region <- 0
    for(j in 1:length(t@x)) {
      if(t@region[j] > originating.region)
        break
      if(previous.region == originating.region & t@region[j] < previous.region) {
        result[t@region[j]] <- result[t@region[j]] + 1
        total.regressions <- total.regressions + 1
        break
      }
      previous.region <- t@region[j]
    }
  }
  if(type=="proportion")
    return(result/total.regressions)
  else if(type=="count")
    return(result)
  else if(type=="posterior.ci.lower")
    return(qbeta(0.025,result + 1, total.regressions - result + 1))
  else if(type=="posterior.ci.upper")
    return(qbeta(0.975,result + 1, total.regressions - result + 1))
  else
    warning(paste("invalid type", type, " -- nothing returned from regression.target.breakdown."))
}

regression.targets.as.data.frame <- function(expt, originating.region) {
  subj <- c()
  item <- c()
  cond <- c()
  outcomes <- c()
  for(i in 1:length(expt@trials)) {
    t <-  expt@trials[[i]]
    highest.region <- 0
    previous.region <- 0
    for(j in 1:length(t@x)) {
      if(t@region[j] > originating.region)
        break
      if(previous.region == originating.region & t@region[j] < previous.region) {
        subj <- c(subj,t@subj)
        item <- c(item,t@item)
        cond <- c(cond,t@condition)
        outcomes <- c(outcomes, t@region[j])
        break
      }
      previous.region <- t@region[j]
    }
  }
  result <- data.frame(subj=subj,item=item,cond=cond,outcome=outcomes)
  return(standardize.df(result))
}


## only returns those trials for which there was a first-pass regression from
## fromRegion. Returns FALSE for absence of such an unbroken regression; TRUE
## for presence
##
## gotta continue to check correctness, and confirm equality conditions on
## delimiters (14 Dec 2008, RPL)
unbroken.fp.regression.to.region <- function(expt,fromRegion,toRegion,only.for.fp.regressions=FALSE,forward.within.region.ok=FALSE) {
  if(toRegion >= fromRegion)
    warning(paste("error: fromRegion", fromRegion," is not larger than toRegion", toRegion, "!"))
  toRegion <- toRegion+1    ## due to this function being bad right now.
  fromRegion <- fromRegion+1
  subj <- c()
  item <- c()
  cond <- c()
  outcomes <- c()
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    leftEdgeFrom <- ifelse(fromRegion==1, 0,
                           expt@delim[t@item, t@condition, fromRegion-1])
    rightEdgeFrom <- ifelse(fromRegion==(dim(expt@delim)[3] + 1), 999,
                            expt@delim[t@item,t@condition, fromRegion])
    leftEdgeTo <- ifelse(toRegion==1, 0,
                           expt@delim[t@item, t@condition, toRegion-1])
    rightEdgeTo <- ifelse(toRegion==(dim(expt@delim)[3] + 1), 999,
                            expt@delim[t@item,t@condition, toRegion])
    passedLeftEdgeFrom <- t@x >= leftEdgeFrom
    passedRightEdgeFrom <- t@x >= rightEdgeFrom
    inFromRegion <- FALSE
    regressedFromFromRegion <- FALSE
    this.outcome <- FALSE
    last.x <- -1
    last.region <- -1
    for(j in 1:length(t@x)) {
      if(inFromRegion & passedRightEdgeFrom[j]) {
        regressedFromFromRegion <- FALSE
        break
      }
      if(inFromRegion & ! passedLeftEdgeFrom[j]) {
        regressedFromFromRegion <- TRUE
        inFromRegion <- FALSE
      }
      if(t@x[j] > leftEdgeFrom & t@x[j] < rightEdgeFrom) {
        inFromRegion <- TRUE
      }
      if(regressedFromFromRegion & (t@x[j] > last.x | t@x[j] < leftEdgeTo)) {
        this.outcome <- FALSE
        break
      }
      if(regressedFromFromRegion & (t@x[j] >= leftEdgeTo & t@x[j] < rightEdgeTo)) {
        this.outcome <- TRUE
        break
      }
      last.x <- t@x[j]
    }
    if(regressedFromFromRegion | ! only.for.fp.regressions) {
      subj <- c(subj, t@subj)
      item <- c(item, t@item)
      cond <- c(cond, t@condition)
      outcomes <- c(outcomes,this.outcome)
    }
  }
  result <- data.frame(subj=subj,item=item,cond=cond,outcome=outcomes)
  return(standardize.df(result))
}

regressedFromYToXWithoutPassingY <- function(expt,x,y,spacesToLeft=FALSE) {
  if(x >= y)
    warning(paste("error: fromRegion", y," is not larger than toRegion", x, "!"))
  subj <- c()
  item <- c()
  cond <- c()
  ord <- c()
  outcomes <- c()
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    regressedToX <- FALSE
    this.outcome <- FALSE
    max.region <- -1
    for(j in 1:length(t@x)) {
      if(spacesToLeft) {
        if(t@x[j] > expt@delim[t@item,t@condition,y])
          break
        else
          if(max.region==y & t@x[j] > expt@delim[t@item,t@condition,x] & t@x[j] <= expt@delim[t@item,t@condition,x+1]) {
            regressedToX <- TRUE
            break
          }
        max.region <- max(max.region, match(TRUE, t@x[j] <= c(expt@delim[t@item,t@condition,],1000)))
      }
      else {
        if(t@region[j] > y)
          break
        else
          if(max.region==y & t@region[j] == x & t@x[j]) {
            regressedToX <- TRUE
            break
          }
        max.region <- max(max.region,t@region[j])
      }
    }
    subj <- c(subj,t@subj)
    item <- c(item,t@item)
    cond <- c(cond,t@condition)
    ord <- c(ord,t@order)
    outcomes <- c(outcomes,regressedToX)
  }
  result <- data.frame(subj=subj,item=item,cond=cond,order=ord,outcome=outcomes)
  return(standardize.df(result))
}

regressedFromYToXWithoutPassingYTime <- function(expt,x,y,spacesToLeft=FALSE) {
  if(x >= y)
    warning(paste("error: fromRegion", y," is not larger than toRegion", x, "!"))
  subj <- c()
  item <- c()
  cond <- c()
  ord <- c()
  rt <- c()
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
#    regressedToX <- FALSE
	this.rt <- 0
    max.region <- -1
    for(j in 1:length(t@x)) {
      if(spacesToLeft) {
        if(t@x[j] > expt@delim[t@item,t@condition,y])
          break
        else
          if(max.region==y & t@x[j] > expt@delim[t@item,t@condition,x] & t@x[j] <= expt@delim[t@item,t@condition,x+1]) 			{
            this.rt <- this.rt + t@endTime[j] - t@startTime[j]
          }
        max.region <- max(max.region, match(TRUE, t@x[j] <= c(expt@delim[t@item,t@condition,],1000)))
      }
      else {
        if(t@region[j] > y)
          break
        else
          if(max.region==y & t@region[j] == x & t@x[j]) {
            this.rt <- this.rt + t@endTime[j] - t@startTime[j]
          }
        max.region <- max(max.region,t@region[j])
      }
    }
    subj <- c(subj,t@subj)
    item <- c(item,t@item)
    cond <- c(cond,t@condition)
    ord <- c(ord,t@order)
    rt <- c(rt,this.rt)
  }
  result <- data.frame(subj=subj,item=item,cond=cond,order=ord,rt=rt)
  return(standardize.df(result))
}

maxCharactersRegressedDuringGoPastReading <- function(expt,y,x=1) {
  subj <- c()
  item <- c()
  cond <- c()
  outcomes <- c()
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    this.outcome <- 0
    max.region <- -1
    y.delim <- expt@delim[t@item,t@condition,y]
    x.delim <- expt@delim[t@item,t@condition,x]
    for(j in 1:length(t@x)) {
      if(t@region[j] > y)
        break
      if(max.region==y & t@region[j] >=x  ) {
        this.outcome <- max(this.outcome,y.delim - t@x[j])
      }
      max.region <- max(max.region,t@region[j])
    }
    this.outcome <- min(this.outcome,y.delim-x.delim)
    subj <- c(subj,t@subj)
    item <- c(item,t@item)
    cond <- c(cond,t@condition)
    outcomes <- c(outcomes,this.outcome)
  }
  result <- data.frame(subj=subj,item=item,cond=cond,outcome=outcomes)
  return(standardize.df(result))
}

## first passes only
regressedFromY <- function(expt,y) {
  subj <- c()
  item <- c()
  cond <- c()
  outcomes <- c()
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    regressed <- FALSE
    max.region <- -1
    for(j in 1:length(t@x)) {
      if(t@region[j] > y)
        break
      if(max.region==y & t@region[j] < y) {
        regressed <- TRUE
        break
      }
      max.region <- max(max.region,t@region[j])
    }
    subj <- c(subj,t@subj)
    item <- c(item,t@item)
    cond <- c(cond,t@condition)
    outcomes <- c(outcomes,regressed)
  }
  result <- data.frame(subj=subj,item=item,cond=cond,outcome=outcomes)
  return(standardize.df(result))
}

## first passes only
regressedFromYtoX <- function(expt,x,y) {
  subj <- c()
  item <- c()
  cond <- c()
  outcomes <- c()
  if(max(x) >= y)
    warning("Can't have x greater than y!")
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    regressed <- NA
    max.region <- -1
    for(j in 1:length(t@x)) {
      if(t@region[j] > y) {
        regressed <- FALSE
        break
      }
      if(max.region==y & t@region[j] < y) {
        if(t@region[j] %in% x)
          regressed <- TRUE
        else
          regressed <- FALSE
        break
      }
      max.region <- max(max.region,t@region[j])
    }
    subj <- c(subj,t@subj)
    item <- c(item,t@item)
    cond <- c(cond,t@condition)
    outcomes <- c(outcomes,regressed)
  }
  result <- data.frame(subj=subj,item=item,cond=cond,outcome=outcomes)
  return(standardize.df(result))
}

## y  is region of interest
## if onlyAcrossRegion is TRUE, then only regressive saccades crossing region count.
## if onlyAcrossRegion is FALSE, then within-region regressive saccades count too.
regressiveSaccadeSize <- function(expt,y, onlyAcrossRegion=TRUE) {
  subj <- c()
  item <- c()
  cond <- c()
  outcomes <- c()
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    outcome <- NA
    max.region <- -1
    previous.region <- -1
    for(j in 1:length(t@x)) {
      if(t@region[j] > y)
        break
      if(j > 1 && previous.region == y && ((onlyAcrossRegion && t@region[j] < y) || ( (! onlyAcrossRegion) && t@x[j] < t@x[j-1]))) {
        outcome <- t@x[j-1] - t@x[j]
        break
      }
      max.region <- max(max.region,t@region[j])
      previous.region <- t@region[j]
    }
    subj <- c(subj,t@subj)
    item <- c(item,t@item)
    cond <- c(cond,t@condition)
    outcomes <- c(outcomes,outcome)
  }
  result <- data.frame(subj=subj,item=item,cond=cond,outcome=outcomes)
  return(standardize.df(result))
}

gotToRegionYWithoutSeeingRegionX <- function(expt,x,y) {
  subj <- c()
  item <- c()
  cond <- c()
  outcomes <- c()
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    outcome <- NA
    for(j in 1:length(t@x)) {
      if(t@region[j] >= y) {
        outcome <- TRUE
        break
      }
      if(t@region[j] ==x) {
        outcome <-  FALSE
        break
      }
    }
    subj <- c(subj,t@subj)
    item <- c(item,t@item)
    cond <- c(cond,t@condition)
    outcomes <- c(outcomes,outcome)
  }
  result <- data.frame(subj=subj,item=item,cond=cond,outcome=outcomes)
  return(standardize.df(result))
}

getTrials <- function(expt,subj,item) {
  result <- list()
  for(i in 1:length(expt@trials)) {
    t <- expt@trials[[i]]
    if(t@subj==subj & t@item==item)
      result <- append(result,t)
  }
  return(result)
}

allSaccadeMatrix <- function(expt,ncond) {
  n <- dim(expt@delim)[3]
  result <- list()
  for(i in 1:ncond) 
    result[[i]] <- matrix(rep(0,n*n),n,n)
  for(t in expt@trials) {
    from <- t@region[-length(t@region)]
    to <- t@region[-1]
    ft <- cbind(from,to)
    result[[t@condition]][ft] <- result[[t@condition]][ft] + 1
  }
  return(result)
}
