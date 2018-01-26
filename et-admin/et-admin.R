### get data to analyze blinks and display change problems in an easy exclusion format

determine.list <- function(df, num.conditions) {
    ## this function assumes that df contains one column for item and
    ## one for condition (which may be NA). It returns a list number,
    ## which is simply the condition number for items 1,
    ## 1+num.conditions, 1+2*num.conditions, ...
    max.item <- max(df$item)
    df <- subset(df, item %in% seq(1, max.item+1, by=num.conditions))
    return(mean(df$cond, na.rm=T))
}

get.correct.cond <- function(item, list, num.conditions) {
    ## Returns the condition for a particular item in a particular
    ## list (assuming standard item rotation)
    return(1 + (item-1+list-1) %% num.conditions)
}

get.timing.exclusions <- function(timing.filename, # output of jhook
                                  timing.max, # maximum ms post fixation onset
                                  asclist.filename, # input to jhook
                                  path="./") { # interpret filenames relative to this
    ## Determines which trials to exclude for various display change
    ## problems, most importantly including having a jhook and
    ## completing too late into fixation. Also excludes for other rare
    ## problems that may indicate display change errors, including
    ## having negative timing values, NA timing values, negative preX
    ## values, having two lines of information produced by jhook, and
    ## completing far too early.

    df.jhook <- read.delim(paste0(path, timing.filename))
    df.jhook$id <- with(df.jhook, paste(subject, item))
    df.jhook$item <- factor(df.jhook$item)
    df.jhook$subject <- factor(df.jhook$subject)
    df.jhook$condition <- factor(df.jhook$condition)
    num.subjs <- length(levels(df.jhook$subject))
    num.items <- length(levels(df.jhook$item))

    ## first, determine what has multiple rows
    df.jhook$dummy <- 1
    df.rows <- with(df.jhook, aggregate(list(count=dummy), list(subject=subject, item=item), sum))
    df.rows$justone <- df.rows$count == 1
    df.rows$id <- with(df.rows, paste(subject, item))
    df.jhook$justone <- df.rows$justone[match(df.jhook$id, df.rows$id)]

    ## now remove everything
    df.jhook <- subset(df.jhook, justone==T) # multiple rows in jhook
                                             # output
    df.jhook <- subset(df.jhook, preX > 0) # impossible
    df.jhook <- subset(df.jhook, !is.na(timing))
    df.jhook <- subset(df.jhook, timing > -40) # shouldn't be able to
                                               # be triggered and
                                               # complete so early...
    df.jhook <- subset(df.jhook, timing <= timing.max)
    df.jhook <- subset(df.jhook, postX < 0) # jhooks
    df.jhook <- subset(df.jhook, ! (comment %in% c("blink", "blink2")))

    ## prepare new dataframe of subj-item combinations
    subj.labels <- read.table(paste0(path, asclist.filename), header=F, col.names=c("label"))
    stopifnot(nrow(subj.labels)==num.subjs)
    subj.labels$number <- seq_len(num.subjs)
    subj.labels$label <- sub("\\", "/", subj.labels$label, fixed=T)
    subj.labels$label <- paste0(path, subj.labels$label)

    df.dispchange <- expand.grid(list(item=seq_len(num.items), subj=seq_len(num.subjs)))
    df.dispchange$id <- with(df.dispchange, paste(subj, item))
    df.dispchange$label <- subj.labels$label[match(df.dispchange$subj, subj.labels$number)]
    df.dispchange$good <- df.dispchange$id %in% df.jhook$id
    df.dispchange$timing <- df.jhook$timing[match(df.dispchange$id, df.jhook$id)]

    ## cleanup and print
    df.jhook$dummy <- NULL
    df.jhook$id <- NULL
    df.jhook$justone <- NULL
    cat("\nHere's the jhook data after exclusions. Please inspect!\n\n")
    print(summary(df.jhook))

    return(df.dispchange[c("subj", "item", "label", "id", "good", "timing")])
}

get.blink.data <- function(da1.filename, num.items) {
    num.cols <- max(count.fields(da1.filename, sep='\t'))
    #df <- read.delim(da1.filename, header=F, col.names=c('trial','cond','item',4:num.cols))
    df <- read.delim(da1.filename, header=F, col.names=seq_len(num.cols))
    df <- subset(df, select=c("X2", "X3"))
    names(df) <- c("cond", "item")
    num.conditions <- max(df$cond)
    list.num <- determine.list(df, num.conditions)
    df$good <- TRUE
    for (item in seq_len(num.items)) {
        correct.cond <- get.correct.cond(item, list.num, num.conditions)
        if (item %in% df$item) {
            # sometimes there are more than one of these (although
            # there shouldn't be and we'll exclude them later)
            stopifnot(all(df$cond[df$item==item]==correct.cond)) 
        } else {
            newrow <- c(correct.cond, item, FALSE)
            df <- rbind(df, newrow)
        }
    }
    df$label <- da1.filename
    return(df)
}

get.da1.list <- function(da1list.filename, path="./") {
    da1s <- read.delim(paste0(path, da1list.filename), header=F, as.is=T)
    da1s <- sub("\\", "/", da1s$V1, fixed=T)
    da1s <- paste0(path, da1s)
}    

get.blink.exclusions <- function(da1list.filename,
                                 num.items,
                                 path="./") {
    subjs <- get.da1.list(da1list.filename, path)
    subj.nums <- seq_along(subjs)
    blinks <- lapply(subjs, function (x) get.blink.data(x, num.items))
    blinks <- do.call(rbind, blinks)
    blinks$subj <- subj.nums[match(blinks$label, subjs)]
    blinks$id <- with(blinks, paste(subj, item))

    blinks$dummy <- 1
    df.rows <- with(blinks, aggregate(list(count=dummy), list(subj=subj, item=item), sum))
    df.rows$justone <- df.rows$count == 1
    df.rows$id <- with(df.rows, paste(subj, item))
    blinks$justone <- df.rows$justone[match(blinks$id, df.rows$id)]
    blinks$good[blinks$justone == F] <- F
    blinks <- subset(blinks, !duplicated(blinks))

    ## cleanup
    blinks$justone <- NULL
    blinks$dummy <- NULL
    blinks <- blinks[with(blinks, order(subj, item)), ]
    row.names(blinks) <- seq_len(nrow(blinks))
    blinks$item <- as.integer(blinks$item)
    return(blinks[c("subj", "item", "label", "id", "good", "cond")])
}

get.expt.from.da1list <- function(da1list.filename,
                                  count.filename,
                                  path="./",
                                  min.fix=80,
                                  ignore.num.fixes=F,
                                  mc.cores=getOption("mc.cores", 2L)) {
    da1s <- get.da1.list(da1list.filename, path)
    subj.nums <- seq_along(da1s)
    expt <- getETexpt(da1s, subj.nums,
                      paste0(path, count.filename),
                      min.fix=min.fix,
                      ignore.num.fixes=ignore.num.fixes,
                      mc.cores=mc.cores)
    return(expt)
}

read.subj.by.item <- function(filename,
                             path="./") {
    df <- read.table(paste0(path, filename))
    num.regions <- ncol(df)-3
    regions <- paste0('region', seq_len(num.regions))
    names(df) <- c('subj', 'item', 'cond', regions)
    return(df)
}
