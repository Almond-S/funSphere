
# copy of mgcv:gam, only with absorb.cons = FALSE passed to gam.setup
#' @import mgcv
gam_ <- function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action,offset=NULL,
                method="GCV.Cp",optimizer=c("outer","newton"),control=list(),#gam.control(),
                scale=0,select=FALSE,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,fit=TRUE,
                paraPen=NULL,G=NULL,in.out=NULL,drop.unused.levels=TRUE,drop.intercept=NULL,
                nei=NULL,discrete=FALSE,...) {
  # attach(getNamespace("mgcv"))
  ## Routine to fit a GAM to some data. The model is stated in the formula, which is then
  ## interpreted to figure out which bits relate to smooth terms and which to parametric terms.
  ## Basic steps:
  ## 1. Formula is split up into parametric and non-parametric parts,
  ##    and a fake formula constructed to be used to pick up data for
  ##    model frame. pterms "terms" object(s) created for parametric
  ##    components, model frame created along with terms object.
  ## 2. 'gam.setup' called to do most of basis construction and other
  ##    elements of model setup.
  ## 3. 'estimate.gam' is called to estimate the model. This performs further
  ##    pre- and post- fitting steps and calls either 'gam.fit' (performance
  ##    iteration) or 'gam.outer' (default method). 'gam.outer' calls the actual
  ##    smoothing parameter optimizer ('newton' by default) and then any post
  ##    processing. The optimizer calls 'gam.fit3/4/5' to estimate the model
  ##    coefficients and obtain derivatives w.r.t. the smoothing parameters.
  ## 4. Finished 'gam' object assembled.
  control <- do.call("gam.control",control)
  if (is.null(G) && discrete) { ## get bam to do the setup
    cl <- match.call() ## NOTE: check all arguments more carefully
    cl[[1]] <- quote(bam)
    cl$fit = FALSE
    G <- eval(cl,parent.frame()) ## NOTE: cl probaby needs modifying in G to work properly (with fit=FALSE reset?? also below??)
  }
  if (is.null(G)) {
    ## create model frame.....
    gp <- interpret.gam(formula) # interpret the formula
    cl <- match.call() # call needed in gam object for update to work
    mf <- match.call(expand.dots=FALSE)
    mf$formula <- gp$fake.formula
    mf$family <- mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp<-mf$H<-mf$select <- mf$drop.intercept <- mf$nei <-
      mf$gamma<-mf$method<-mf$fit<-mf$paraPen<-mf$G<-mf$optimizer <- mf$in.out <- mf$discrete <- mf$...<-NULL
    mf$drop.unused.levels <- drop.unused.levels
    mf[[1]] <- quote(stats::model.frame) ## as.name("model.frame")
    pmf <- mf
    mf <- eval(mf, parent.frame()) # the model frame now contains all the data
    if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
    terms <- attr(mf,"terms")

    ## summarize the *raw* input variables
    ## note can't use get_all_vars here -- buggy with matrices
    vars <- all.vars1(gp$fake.formula[-2]) ## drop response here
    inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))

    ## allow a bit of extra flexibility in what `data' is allowed to be (as model.frame actually does)
    if (!is.list(data)&&!is.data.frame(data)) data <- as.data.frame(data)

    dl <- eval(inp, data, parent.frame())
    names(dl) <- vars ## list of all variables needed
    var.summary <- variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
    rm(dl) ## save space

    ## pterms are terms objects for the parametric model components used in
    ## model setup - don't try obtaining by evaluating pf in mf - doesn't
    ## work in general (e.g. with offset)...

    if (is.list(formula)) { ## then there are several linear predictors
      environment(formula) <- environment(formula[[1]]) ## e.g. termplots needs this
      pterms <- list()
      tlab <- rep("",0)
      for (i in 1:length(formula)) {
        pmf$formula <- gp[[i]]$pf
        pterms[[i]] <- attr(eval(pmf, parent.frame()),"terms")
        tlabi <- attr(pterms[[i]],"term.labels")
        if (i>1&&length(tlabi)>0) tlabi <- paste(tlabi,i-1,sep=".")
        tlab <- c(tlab,tlabi)
      }
      attr(pterms,"term.labels") <- tlab ## labels for all parametric terms, distinguished by predictor
    } else { ## single linear predictor case
      pmf$formula <- gp$pf
      pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part
      pterms <- attr(pmf,"terms") ## pmf only used for this
    }

    if (is.character(family)) family <- eval(parse(text=family))
    if (is.function(family)) family <- family()
    if (is.null(family$family)) stop("family not recognized")

    if (family$family[1]=="gaussian" && family$link=="identity") am <- TRUE
    else am <- FALSE

    if (!control$keepData) rm(data) ## save space

    ## check whether family requires intercept to be dropped...
    # drop.intercept <- if (is.null(family$drop.intercept) || !family$drop.intercept) FALSE else TRUE
    # drop.intercept <- as.logical(family$drop.intercept)
    if (is.null(family$drop.intercept)) { ## family does not provide information
      lengthf <- if (is.list(formula)) length(formula) else 1
      if (is.null(drop.intercept)) drop.intercept <- rep(FALSE, lengthf) else {
        drop.intercept <- rep(drop.intercept,length=lengthf) ## force drop.intercept to correct length
        if (sum(drop.intercept)) family$drop.intercept <- drop.intercept ## ensure prediction works
      }
    } else drop.intercept <- as.logical(family$drop.intercept) ## family overrides argument

    if (inherits(family,"general.family")&&!is.null(family$presetup)) eval(family$presetup)

    gsname <- if (is.list(formula)) "gam.setup.list" else "gam.setup"

    G <- do.call(gsname,list(formula=gp,pterms=pterms,
                             data=mf,knots=knots,sp=sp,min.sp=min.sp,
                             H=H,absorb.cons=FALSE,sparse.cons=0,select=select,
                             idLinksBases=control$idLinksBases,scale.penalty=control$scalePenalty,
                             paraPen=paraPen,drop.intercept=drop.intercept))

    G$var.summary <- var.summary
    G$family <- family

    if ((is.list(formula)&&(is.null(family$nlp)||family$nlp!=gp$nlp))||
        (!is.list(formula)&&!is.null(family$npl)&&(family$npl>1))) stop("incorrect number of linear predictors for family")

    G$terms<-terms;
    G$mf<-mf;G$cl<-cl;
    G$am <- am

    if (is.null(G$offset)) G$offset<-rep(0,G$n)

    G$min.edf <- G$nsdf ## -dim(G$C)[1]
    if (G$m) for (i in 1:G$m) G$min.edf<-G$min.edf+G$smooth[[i]]$null.space.dim

    G$formula <- formula
    G$pred.formula <- gp$pred.formula
    environment(G$formula)<-environment(formula)
  } else { ## G not null
    if (!is.null(sp)&&any(sp>=0)) { ## request to modify smoothing parameters
      if (is.null(G$L)) G$L <- diag(length(G$sp))
      if (length(sp)!=ncol(G$L)) stop('length of sp must be number of free smoothing parameters in original model')
      ind <- sp>=0 ## which smoothing parameters are now fixed
      spind <- log(sp[ind]);
      spind[!is.finite(spind)] <- -30 ## set any zero parameters to effective zero
      G$lsp0 <- G$lsp0 + drop(G$L[,ind,drop=FALSE] %*% spind) ## add fix to lsp0
      G$L <- G$L[,!ind,drop=FALSE] ## drop the cols of G
      G$sp <- rep(-1,ncol(G$L))
    }
  }

  if (!fit) {
    class(G) <- "gam.prefit"
    return(G)
  }

  if (ncol(G$X)>nrow(G$X)) stop("Model has more coefficients than data")

  G$conv.tol <- control$mgcv.tol      # tolerence for mgcv
  G$max.half <- control$mgcv.half # max step halving in Newton update mgcv

  object <- estimate.gam(G,method,optimizer,control,in.out,scale,gamma,nei=nei,...)


  if (!is.null(G$L)) {
    object$full.sp <- as.numeric(exp(G$L%*%log(object$sp)+G$lsp0))
    names(object$full.sp) <- names(G$lsp0)
  }
  names(object$sp) <- names(G$sp)
  object$paraPen <- G$pP
  object$formula <- G$formula
  ## store any lpi attribute of G$X for use in predict.gam...
  if (is.list(object$formula)) attr(object$formula,"lpi") <- attr(G$X,"lpi")
  object$var.summary <- G$var.summary
  object$cmX <- G$cmX ## column means of model matrix --- useful for CIs
  object$model<-G$mf # store the model frame
  object$na.action <- attr(G$mf,"na.action") # how to deal with NA's
  object$control <- control
  object$terms <- G$terms
  object$pred.formula <- G$pred.formula

  attr(object$pred.formula,"full") <- reformulate(all.vars(object$terms))

  object$pterms <- G$pterms
  object$assign <- G$assign # applies only to pterms
  object$contrasts <- G$contrasts
  object$xlevels <- G$xlevels
  object$offset <- G$offset
  if (!is.null(G$Xcentre)) object$Xcentre <- G$Xcentre
  if (control$keepData) object$data <- data
  object$df.residual <- nrow(G$X) - sum(object$edf)
  object$min.edf <- G$min.edf
  if (G$am&&!(method%in%c("REML","ML","P-ML","P-REML"))) object$optimizer <- "magic" else object$optimizer <- optimizer
  object$call <- G$cl # needed for update() to work
  class(object) <- c("gam","glm","lm")
  if (is.null(object$deviance)) object$deviance <- sum(residuals(object,"deviance")^2)
  names(object$gcv.ubre) <- method
  ## The following lines avoid potentially very large objects in hidden environments being stored
  ## with fitted gam objects. The downside is that functions like 'termplot' that rely on searching in
  ## the environment of the formula can fail...
  environment(object$formula) <- environment(object$pred.formula) <-
    environment(object$terms) <- environment(object$pterms) <- .GlobalEnv
  if (!is.null(object$model))  environment(attr(object$model,"terms"))  <- .GlobalEnv
  if (!is.null(attr(object$pred.formula,"full"))) environment(attr(object$pred.formula,"full")) <- .GlobalEnv
  object
} ## gam


# copy of mgcv:bam, only with absorb.cons = FALSE passed to gam.setup
#' @import mgcv
bam_ <- function(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,na.action=na.omit,
                offset=NULL,method="fREML",control=list(),select=FALSE,scale=0,gamma=1,knots=NULL,sp=NULL,
                min.sp=NULL,paraPen=NULL,chunk.size=10000,rho=0,AR.start=NULL,discrete=FALSE,
                cluster=NULL,nthreads=1,gc.level=0,use.chol=FALSE,samfrac=1,coef=NULL,
                drop.unused.levels=TRUE,G=NULL,fit=TRUE,drop.intercept=NULL,...)
  # attach(getNamespace("mgcv"))
  ## Routine to fit an additive model to a large dataset. The model is stated in the formula,
  ## which is then interpreted to figure out which bits relate to smooth terms and which to
  ## parametric terms.
  ## This is a modification of `gam' designed to build the QR decomposition of the model matrix
  ## up in chunks, to keep memory costs down.
  ## If cluster is a parallel package cluster uses parallel QR build on cluster.
  ## 'n.threads' is number of threads to use for non-cluster computation (e.g. combining
  ## results from cluster nodes). If 'NA' then is set to max(1,length(cluster)).
{ control <- do.call("gam.control",control)
if (control$trace) t3 <- t2 <- t1 <- t0 <- proc.time()
if (length(nthreads)==1) nthreads <- rep(nthreads,2)
if (is.null(G)) { ## need to set up model!
  if (is.character(family))
    family <- eval(parse(text = family))
  if (is.function(family))
    family <- family()
  if (is.null(family$family))
    stop("family not recognized")

  if (family$family=="gaussian"&&family$link=="identity") am <- TRUE else am <- FALSE
  if (scale==0) { if (family$family%in%c("poisson","binomial")) scale <- 1 else scale <- -1}
  if (!method%in%c("fREML","GACV.Cp","GCV.Cp","REML",
                   "ML","P-REML","P-ML")) stop("un-supported smoothness selection method")
  if (is.logical(discrete)) {
    discretize <- discrete
    discrete <- NULL ## use default discretization, if any
  } else {
    discretize <- if (is.numeric(discrete)) TRUE else FALSE
  }
  if (discretize) {
    if (method!="fREML") {
      discretize <- FALSE
      warning("discretization only available with fREML")
    } else {
      if (!is.null(cluster)) warning("discrete method does not use parallel cluster - use nthreads instead")
      if (all(is.finite(nthreads)) && any(nthreads>1) && !mgcv.omp()) warning("openMP not available: single threaded computation only")
    }
  }
  if (inherits(family,"extended.family")) {
    family <- fix.family.link(family); efam <- TRUE
  } else efam <- FALSE

  if (method%in%c("fREML")&&!is.null(min.sp)) {
    min.sp <- NULL
    warning("min.sp not supported with fast REML computation, and ignored.")
  }

  gp <- interpret.gam(formula) # interpret the formula
  if (discretize && length(gp$smooth.spec)==0) {
    ok <- TRUE
    ## check it's not a list formula
    if (!is.null(gp$nlp)) for (i in 1:gp$nlp) if (length(gp[[i]]$smooth.spec)>0) ok <- FALSE
    if (ok) {
      warning("no smooths, ignoring `discrete=TRUE'")
      discretize <- FALSE
    }
  }
  if (discretize) {
    ## re-order the tensor terms for maximum efficiency, and
    ## signal that "re"/"fs" terms should be constructed with marginals
    ## also for efficiency

    if (is.null(gp$nlp)) for (i in 1:length(gp$smooth.spec)) {
      if (inherits(gp$smooth.spec[[i]],"tensor.smooth.spec")) gp$smooth.spec[[i]] <- tero(gp$smooth.spec[[i]])
      #if (inherits(gp$smooth.spec[[i]],c("re.smooth.spec","fs.smooth.spec"))&&gp$smooth.spec[[i]]$dim>1)
      if (!is.null(gp$smooth.spec[[i]]$tensor.possible)&&gp$smooth.spec[[i]]$dim>1){
        class(gp$smooth.spec[[i]]) <- c(class(gp$smooth.spec[[i]]),"tensor.smooth.spec")
        if (is.null(gp$smooth.spec[[i]]$margin)) {
          gp$smooth.spec[[i]]$margin <- list()
          ## only ok for 'fs' with univariate metric variable (caught in 'fs' construcor)...
          for (j in 1:gp$smooth.spec[[i]]$dim) gp$smooth.spec[[i]]$margin[[j]] <- list(term=gp$smooth.spec[[i]]$term[j])
        }
      }
    } else for (j in 1:length(formula)) if (length(gp[[j]]$smooth.spec)>0) for (i in 1:length(gp[[j]]$smooth.spec)) {
      if (inherits(gp[[j]]$smooth.spec[[i]],"tensor.smooth.spec")) gp[[j]]$smooth.spec[[i]] <- tero(gp[[j]]$smooth.spec[[i]])
      #if (inherits(gp[[j]]$smooth.spec[[i]],c("re.smooth.spec","fs.smooth.spec"))&&gp[[j]]$smooth.spec[[i]]$dim>1)
      if (!is.null(gp[[j]]$smooth.spec[[i]]$tensor.possible)&&gp[[j]]$smooth.spec[[i]]$dim>1) {
        class(gp[[j]]$smooth.spec[[i]]) <- c(class(gp[[j]]$smooth.spec[[i]]),"tensor.smooth.spec")
        if (is.null(gp[[j]]$smooth.spec[[i]]$margin)) {
          gp[[j]]$smooth.spec[[i]]$margin <- list()
          ## only ok for 'fs' with univariate metric variable (caught in 'fs' construcor)...
          for (k in 1:gp[[j]]$smooth.spec[[i]]$dim) gp[[j]]$smooth.spec[[i]]$margin[[k]] <- list(term=gp[[j]]$smooth.spec[[i]]$term[k])
        }
      }
    }
  } ## if (discretize)
  cl <- match.call() # call needed in gam object for update to work
  mf <- match.call(expand.dots=FALSE)
  mf$formula <- gp$fake.formula
  mf$method <-  mf$family<-mf$control<-mf$scale<-mf$knots<-mf$sp<-mf$min.sp <- mf$gc.level <-
    mf$gamma <- mf$paraPen<- mf$chunk.size <- mf$rho  <- mf$cluster <- mf$discrete <-
    mf$use.chol <- mf$samfrac <- mf$nthreads <- mf$G <- mf$fit <- mf$select <- mf$drop.intercept <-
    mf$coef <- mf$...<-NULL
  mf$drop.unused.levels <- drop.unused.levels
  mf[[1]] <- quote(stats::model.frame) ## as.name("model.frame")

  if (is.list(formula)) { ## then there are several linear predictors
    environment(formula) <- environment(formula[[1]]) ## e.g. termplots needs this
    pterms <- list()
    tlab <- rep("",0)
    pmf.names <- rep("",0)
    for (i in 1:length(formula)) {
      pmf <- mf
      pmf$formula <- gp[[i]]$pf
      pmf <- eval(pmf, parent.frame())
      pmf.names <- c(pmf.names,names(pmf))
      pterms[[i]] <- attr(pmf,"terms")
      tlabi <- attr(pterms[[i]],"term.labels")
      if (i>1&&length(tlabi)>0) tlabi <- paste(tlabi,i-1,sep=".")
      tlab <- c(tlab,tlabi)
    }
    pmf.names <- unique(pmf.names)
    attr(pterms,"term.labels") <- tlab ## labels for all parametric terms, distinguished by predictor
    nlp <- gp$nlp
    lpid <- list() ## list of terms for each lp
    lpid[[nlp]] <- rep(0,0)
  } else { ## single linear predictor case
    nlp <- 1
    pmf <- mf
    pmf$formula <- gp$pf
    pmf <- eval(pmf, parent.frame()) # pmf contains all data for parametric part
    pterms <- attr(pmf,"terms") ## pmf only used for this and discretization, if selected.
    pmf.names <- names(pmf)
  }

  if (gc.level>0) gc()

  mf <- eval(mf, parent.frame()) # the model frame now contains all the data

  if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
  terms <- attr(mf,"terms")
  if (gc.level>0) gc()
  if (rho!=0&&!is.null(mf$"(AR.start)")) if (!is.logical(mf$"(AR.start)")) stop("AR.start must be logical")

  ## summarize the *raw* input variables
  ## note can't use get_all_vars here -- buggy with matrices
  vars <- all.vars1(gp$fake.formula[-2]) ## drop response here
  inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))

  ## allow a bit of extra flexibility in what `data' is allowed to be (as model.frame actually does)
  if (!is.list(data)&&!is.data.frame(data)) data <- as.data.frame(data)

  dl <- eval(inp, data, parent.frame())
  if (!control$keepData) { rm(data);if (gc.level>0) gc()} ## save space
  names(dl) <- vars ## list of all variables needed
  var.summary <- variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
  rm(dl); if (gc.level>0) gc() ## save space

  ## should we force the intercept to be dropped, meaning that the constant is removed
  ## from the span of the parametric effects?
  if (is.null(family$drop.intercept)) { ## family does not provide information
    lengthf <- if (is.list(formula)) length(formula) else 1
    if (is.null(drop.intercept)) drop.intercept <- rep(FALSE,lengthf) else {
      drop.intercept <- rep(drop.intercept,length=lengthf) ## force drop.intercept to correct length
      if (sum(drop.intercept)) family$drop.intercept <- drop.intercept ## ensure prediction works
    }
  } else drop.intercept <- as.logical(family$drop.intercept) ## family overrides argument


  ## need mini.mf for basis setup, then accumulate full X, y, w and offset
  if (discretize) {
    ## discretize the data, creating list mf0 with discrete values
    ## and indices giving the discretized value for each element of model frame.
    ## 'discrete' can be null, or contain a discretization size, or
    ## a discretization size per smooth term.
    dk <- discrete.mf(gp,mf,pmf.names,m=discrete)
    mf0 <- dk$mf ## padded discretized model frame
    sparse.cons <- 0 ## default constraints required for tensor terms

  } else {
    mf0 <- mini.mf(mf,chunk.size)
    sparse.cons <- -1
  }
  rm(pmf); ## no further use

  ## allow bam to set up general families, even if it can not (yet) process them...
  if (inherits(family,"general.family")&&!is.null(family$presetup)) eval(family$presetup)
  gsname <- if (is.list(formula)) "gam.setup.list" else "gam.setup"

  if (control$trace) t1 <- proc.time()
  reset <- TRUE
  while (reset) {
    G <- do.call(gsname,list(formula=gp,pterms=pterms,
                             data=mf0,knots=knots,sp=sp,min.sp=min.sp,
                             H=NULL,absorb.cons=FALSE,sparse.cons=sparse.cons,select=select,
                             idLinksBases=!discretize,scale.penalty=control$scalePenalty,
                             paraPen=paraPen,apply.by=!discretize,drop.intercept=drop.intercept,modCon=2))

    if (!discretize&&ncol(G$X)>=chunk.size) { ## no point having chunk.size < p
      chunk.size <- 4*ncol(G$X)
      warning(gettextf("chunk.size < number of coefficients. Reset to %d",chunk.size))
      if (chunk.size>=nrow(mf)) { ## no sense splitting up computation
        mf0 <- mf ## just use full dataset
      } else reset <- FALSE
    } else reset <- FALSE
  }
  if (control$trace) t2 <- proc.time()
  if (discretize) {
    ks <- matrix(0,0,2) ## NOTE: slightly more efficient not to repeatedly extend
    if (nlp>1) lpi <- attr(G$X,"lpi")
    v <- G$Xd <- list()
    kb <- k <- 1 ## kb indexes blocks, k indexes marginal matrices
    G$kd <- dk$k
    qc <- dt <- ts <- rep(0,length(G$smooth))
    ## have to extract full parametric model matrix from pterms and mf
    npt <- if (nlp==1) 1 else length(G$pterms)
    lpip <- list() ## record coef indices for each discretized term
    for (j in 1:npt) { ## loop over parametric terms in each formula
      paratens <- TRUE
      ## get the parametric model split into tensor components...
      ptens <- if (nlp==1&&!is.list(G$pterms)) terms2tensor(G$pterms,mf0,drop.intercept=drop.intercept[j]) else
        terms2tensor(G$pterms[[j]],mf0,drop.intercept=drop.intercept[j])
      if (j>1) ptens$term.labels <- paste(ptens$term.labels,".",j-1,sep="")
      ## now locate the index vectors for each parametric marginal discrete matrix ptens$X
      if (!is.null(ptens)) {
        np <-length(ptens$X);n <- nrow(mf)
        qc <- c(rep(0,length(ptens$ts)),qc) ## extend (empty) constraint indicator
        coef.ind <-1
        kk <- 1
        for (i in 1:length(ptens$ts)) {
          jj <- 1
          ts[kb] = k;names(ts)[kb] <- ptens$term.labels[i]
          dt[kb] = ptens$dt[i]
          for (ii in 1:dt[kb]) {
            ks <- rbind(ks,dk$ks[ptens$xname[kk],])
            G$Xd[[k]] <- ptens$X[[kk]][1:dk$nr[ptens$xname[kk]],,drop=FALSE];
            kk <- kk + 1
            jj <- jj * ncol(G$Xd[[k]]) ## number of coeffs for this term
            k <- k + 1 ## update matrix counter
          }
          lpip[[kb]] <- 1:jj - 1 + if (nlp==1) coef.ind else attr(G$nsdf,"pstart")[j] - 1 + coef.ind
          coef.ind <- coef.ind + jj
          if (nlp>1) for (ii in 1:length(lpi)) if (any(lpip[[kb]]%in%lpi[[ii]])) lpid[[ii]] <- c(lpid[[ii]],kb)
          kb <- kb + 1 ## update block counter
        }
      } ## is.null(ptens)
    } ## loop over parametric terms in each formula
    ## k is marginal counter, kb is block counter
    ## G$kd[,ks[j,1]:ks[j,2]] (dk$k) gives index columns for term j, thereby allowing
    ## summation over matrix covariates....
    #G$ks <- cbind(dk$k.start[-length(dk$k.start)],dk$k.start[-1])


    drop <- rep(0,0) ## index of te related columns to drop
    if (length(G$smooth)>0) for (i in 1:length(G$smooth)) { ## loop over smooths
      ## potentially each smoother model matrix can be made up of a sequence
      ## of row-tensor products, nead to loop over such sub blocks...
      nsub <- if (!is.null(G$smooth[[i]]$ts)) length(G$smooth[[i]]$ts) else 1
      lp0 <- G$smooth[[i]]$first.para -1  ## offset for start of coeffs for this sub block
      for (sb in 1:nsub) { ## loop over sub-blocks
        np <- 1 ## compute number of sub-block coeffs
        ts[kb] <- k;names(ts)[kb] <- G$smooth[[i]]$label
        ## first deal with any by variable (as first marginal of tensor)...
        if (G$smooth[[i]]$by!="NA") {
          dt[kb] <- 1
          termk <- G$smooth[[i]]$by
          by.var <- dk$mf[[termk]][1:dk$nr[termk]]
          if (is.factor(by.var)) {
            ## create dummy by variable...
            by.var <- as.numeric(by.var==G$smooth[[i]]$by.level)
          }
          G$Xd[[k]] <- matrix(by.var,dk$nr[termk],1)
          np <- ncol(G$Xd[[k]])
          ks <- rbind(ks,dk$ks[termk,])
          k <- k + 1
          by.present <- 1
        } else by.present <- dt[kb] <- 0
        ## ... by done
        if (inherits(G$smooth[[i]],"tensor.smooth")) {
          nmar <- if (is.null(G$smooth[[i]]$dt)) length(G$smooth[[i]]$margin) else G$smooth[[i]]$dt[sb]
          dt[kb] <- dt[kb] + nmar

          jind <- if (sb>1) G$smooth[[i]]$ts[sb] + 1:G$smooth[[i]]$dt[sb] - 1 else 1:nmar
          for (j in jind) {
            termk <- G$smooth[[i]]$margin[[j]]$term[1]
            G$Xd[[k]] <- G$smooth[[i]]$margin[[j]]$X[1:dk$nr[termk],,drop=FALSE]
            np <- np * ncol(G$Xd[[k]])
            ks <- rbind(ks,dk$ks[termk,])
            k <- k + 1
          }
          ## deal with any side constraints on tensor terms
          if (sb==1) { ## only once per smooth!
            di <- attr(G$smooth[[i]],"del.index")
            if (!is.null(di)&&length(di>0)) {
              di <- di + G$smooth[[i]]$first.para + length(drop)  - 1
              drop <- c(drop,di)
            }

            ## deal with tensor smooth constraint
            qrc <- attr(G$smooth[[i]],"qrc")
            ## compute v such that Q = I-vv' and Q[,-1] is constraint null space basis
            if (inherits(qrc,"qr")) {
              v[[kb]] <- qrc$qr/sqrt(qrc$qraux);v[[kb]][1] <- sqrt(qrc$qraux)
              qc[kb] <- 1 ## indicate a constraint
            } else if (length(qrc)>1) { ## Kronecker product of set to zero constraints
              ## on entry qrc is [unused.index, dim1, dim2,..., total number of constraints]
              v[[kb]] <- c(length(qrc)-2,qrc[-1]) ## number of sum-to-zero contrasts, their dimensions, number of constraints
              qc[kb] <- -1
            } else {
              v[[kb]] <- rep(0,0) ##
              if (!inherits(qrc,"character")||qrc!="no constraints") warning("unknown tensor constraint type")
            }
          } else { ## sb==1 once per smooth stuff
            qc <- c(qc,0) ## extend
            v[[kb]] <- rep(0,0)
          }
        } else { ## not a tensor smooth
          v[[kb]] <- rep(0,0)
          dt[kb] <- dt[kb] + 1
          termk <- G$smooth[[i]]$term[1]
          G$Xd[[k]] <- G$X[1:dk$nr[termk],G$smooth[[i]]$first.para:G$smooth[[i]]$last.para,drop=FALSE]
          np <- np * ncol(G$Xd[[k]])
          ks <- rbind(ks,dk$ks[termk,])
          k <- k + 1
        }
        #jj <- G$smooth[[i]]$first.para:G$smooth[[i]]$last.para;
        if (sb==1&&qc[kb]) {
          ncon <- if (qc[kb]<0) v[[kb]][length(v[[kb]])] else 1
          jj <- 1:(np-ncon) + lp0; lp0 <- lp0 + np - ncon
          ## Hard to think of an application requiring constraint when nsub>1, hence not
          ## worked out yet. Add warning to make sure this is flagged if attempt made
          ## to do this in future....
          if (nsub>1) warning("constraints for sub blocked tensor products un-tested")
        } else {
          jj <- 1:np + lp0; lp0 <- lp0 + np
        }
        lpip[[kb]] <- jj
        if (nlp>1) { ## record which lp each discrete term belongs to (can be more than one)
          for (j in 1:nlp) if (any(jj %in% lpi[[j]])) lpid[[j]] <- c(lpid[[j]],kb)
        }
        kb <- kb + 1
      } ## sub block loop
    } ## looping over smooths
    ## put lpid indices into coefficient index order...
    if (nlp>1) {
      for (j in 1:nlp) lpid[[j]] <- lpid[[j]][order(unlist(lapply(lpip[lpid[[j]]],max)))]
      G$lpid <- lpid
    }
    if (length(drop>0)) G$drop <- drop ## index of terms to drop as a result of side cons on tensor terms
    attr(G$Xd,"lpip") <- lpip ## index of coefs by term
    ## ... Xd is the list of discretized model matrices, or marginal model matrices
    ## kd contains indexing vectors, so the ith model matrix or margin is Xd[[i]][kd[i,],]
    ## ts[i] is the starting matrix in Xd for the ith model matrix, while dt[i] is the number
    ## of elements of Xd that make it up (1 for a singleton, more for a tensor).
    ## v is list of Householder vectors encoding constraints and qc the constraint indicator.
    G$v <- v;G$ts <- ts;G$dt <- dt;G$qc <- qc

    G$ks <- ks
    jj <- max(G$ks)-1
    if (ncol(G$kd) > jj) G$kd <- G$kd[,1:jj]
    ## bundle up discretization information needed for discrete prediction...
    G$dinfo <- list(gp=gp, v = G$v, ts = G$ts, dt = G$dt, qc = G$qc, drop = G$drop, pmf.names=pmf.names,lpip=lpip)
    if (nlp>1) G$dinfo$lpid <- lpid
    if (paratens) G$dinfo$para.discrete <- TRUE
  } ## if (discretize)

  if (control$trace) t3 <- proc.time()

  ## no advantage to "fREML" with no free smooths...
  if (((!is.null(G$L)&&ncol(G$L) < 1)||(length(G$sp)==0))&&method=="fREML") method <- "REML"

  G$var.summary <- var.summary
  G$family <- family
  G$terms<-terms;
  G$pred.formula <- gp$pred.formula

  n <- nrow(mf)

  if (is.null(mf$"(weights)")) G$w<-rep(1,n)
  else G$w<-mf$"(weights)"

  G$y <- mf[[gp$response]]
  ## now get offset, dealing with possibility of multiple predictors (see gam.setup)
  ## the point is that G$offset relates to the compressed or discretized model frame,
  ## so we need to correct it to the full data version...
  if (discretize) {
    if (is.list(pterms)) { ## multiple predictors
      for (i in 1:length(pterms)) {
        offi <- attr(pterms[[i]],"offset")
        if (is.null(offi)) G$offset[[i]] <- rep(0,n) else {
          G$offset[[i]] <- mf[[names(attr(pterms[[i]],"dataClasses"))[offi]]]
          if (is.null(G$offset[[i]])) G$offset[[i]] <- rep(0,n)
        }
      }
    } else { ## single predictor, handle as non-discrete
      G$offset <- model.offset(mf)
      if (is.null(G$offset)) G$offset <- rep(0,n)
    }
  } else { ## non-discrete
    G$offset <- model.offset(mf)
    if (is.null(G$offset)) G$offset <- rep(0,n)
  }

  if (!discretize && ncol(G$X)>nrow(mf)) stop("Model has more coefficients than data")

  if (ncol(G$X) > chunk.size && !discretize) { ## no sense having chunk.size < p
    chunk.size <- 4*ncol(G$X)
    warning(gettextf("chunk.size < number of coefficients. Reset to %d",chunk.size))
  }

  G$cl <- cl
  G$am <- am

  G$min.edf<-G$nsdf #-dim(G$C)[1]
  if (G$m) for (i in 1:G$m) G$min.edf<-G$min.edf+G$smooth[[i]]$null.space.dim
  G$discretize <- discretize
  G$formula<-formula
  ## environment(G$formula)<-environment(formula)
  environment(G$pterms) <- environment(G$terms) <- environment(G$pred.formula) <-
    environment(G$formula) <- .BaseNamespaceEnv

} else { ## G supplied
  if (scale<=0) scale <- G$scale
  efam <- G$efam
  mf <- G$mf; G$mf <- NULL
  gp <- G$gp; G$gp <- NULL
  na.action <- G$na.action; G$na.action <- NULL
  if (!is.null(sp)&&any(sp>=0)) { ## request to modify smoothing parameters
    if (is.null(G$L)) G$L <- diag(length(G$sp))
    if (length(sp)!=ncol(G$L)) stop('length of sp must be number of free smoothing parameters in original model')
    ind <- sp>=0 ## which smoothing parameters are now fixed
    spind <- log(sp[ind]);
    spind[!is.finite(spind)] <- -30 ## set any zero parameters to effective zero
    G$lsp0 <- G$lsp0 + drop(G$L[,ind,drop=FALSE] %*% spind) ## add fix to lsp0
    G$L <- G$L[,!ind,drop=FALSE] ## drop the cols of G
    G$sp <- rep(-1,ncol(G$L))
  }
} ## end of G setup

if (!fit) {
  G$efam <- efam
  G$scale <- scale
  G$mf <- mf;G$na.action <- na.action;G$gp <- gp
  class(G) <- "bam.prefit"
  return(G)
}

if (inherits(G$family,"general.family")) stop("general families not supported by bam")

## number of threads to use for non-cluster node computation
if (!is.finite(nthreads[1])||nthreads[1]<1) nthreads[1] <- max(1,length(cluster))

G$conv.tol<-control$mgcv.tol      # tolerence for mgcv
G$max.half<-control$mgcv.half     # max step halving in bfgs optimization


## now build up proper model matrix, and deal with y, w, and offset...

if (control$trace) cat("Setup complete. Calling fit\n")

colnamesX <- colnames(G$X)

if (G$am&&!G$discretize) {
  if (nrow(mf)>chunk.size) G$X <- matrix(0,0,ncol(G$X)); if (gc.level>1) gc()
  object <- bam.fit(G,mf,chunk.size,gp,scale,gamma,method,rho=rho,cl=cluster,
                    gc.level=gc.level,use.chol=use.chol,npt=nthreads[1])
} else if (G$discretize) {
  object <- bgam.fitd(G, mf, gp ,scale ,nobs.extra=0,rho=rho,coef=coef,
                      control = control,npt=nthreads,gc.level=gc.level,gamma=gamma,...)

} else {
  G$X  <- matrix(0,0,ncol(G$X)); if (gc.level>1) gc()
  if (rho!=0) warning("AR1 parameter rho unused with generalized model")
  if (samfrac<1 && samfrac>0) { ## sub-sample first to get close to right answer...
    ind <- sample(1:nrow(mf),ceiling(nrow(mf)*samfrac))
    if (length(ind)<2*ncol(G$X)) warning("samfrac too small - ignored") else {
      Gw <- G$w;Goffset <- G$offset
      G$w <- G$w[ind];G$offset <- G$offset[ind]
      control1 <- control
      control1$epsilon <- 1e-2
      object <- bgam.fit(G, mf[ind,], chunk.size, gp ,scale ,gamma,method=method,nobs.extra=0,
                         control = control1,cl=cluster,npt=nthreads[1],gc.level=gc.level,coef=coef,
                         use.chol=use.chol,samfrac=1,...)
      G$w <- Gw;G$offset <- Goffset
      coef <- object$coefficients
    }
  }
  ## fit full dataset
  object <- bgam.fit(G, mf, chunk.size, gp ,scale ,gamma,method=method,coef=coef,
                     control = control,cl=cluster,npt=nthreads[1],gc.level=gc.level,
                     use.chol=use.chol,...)
}

if (gc.level>0) gc()

if (control$trace) t4 <- proc.time()

if (control$trace) cat("Fit complete. Finishing gam object.\n")

if (scale < 0) { object$scale.estimated <- TRUE;object$scale <- object$scale.est} else {
  object$scale.estimated <- FALSE; object$scale <- scale
}

object$assign <- G$assign # applies only to pterms
object$boundary <- FALSE  # always FALSE for this case
object$call<-G$cl # needed for update() to work
object$cmX <- G$cmX ## column means of model matrix --- useful for CIs

object$contrasts <- G$contrasts
object$control <- control
object$converged <- TRUE ## no iteration
object$data <- NA ## not saving it in this case
object$df.null <- nrow(mf)
object$df.residual <- object$df.null - sum(object$edf)

if (is.null(object$family)) object$family <- family
object$formula <- G$formula

if (method=="GCV.Cp") {
  if (scale<=0) object$method <- "GCV" else object$method <- "UBRE"
} else {
  object$method <- method
}
object$min.edf<-G$min.edf
object$model <- mf;rm(mf);if (gc.level>0) gc()
object$na.action <- attr(object$model,"na.action") # how to deal with NA's
object$nsdf <- G$nsdf
if (G$nsdf>0) names(object$coefficients)[1:G$nsdf] <- colnamesX[1:G$nsdf]
object$offset <- G$offset
##object$prior.weights <- G$w
object$pterms <- G$pterms
object$pred.formula <- G$pred.formula
object$smooth <- G$smooth

object$terms <- G$terms
object$var.summary <- G$var.summary
if (is.null(object$wt)) object$weights <- object$prior.weights else
  object$weights <- object$wt
object$xlevels <- G$xlevels
#object$y <- object$model[[gp$response]]
object$NA.action <- na.action ## version to use in bam.update
names(object$sp) <- names(G$sp)
if (!is.null(object$full.sp)) names(object$full.sp) <- names(G$lsp0)

names(object$coefficients) <- G$term.names
names(object$edf) <- G$term.names

## note that predict.gam assumes that it must be ok not to split the
## model frame, if no new data supplied, so need to supply explicitly
class(object) <- c("bam","gam","glm","lm")
if (!G$discretize) { object$linear.predictors <-
  as.numeric(predict.bam(object,newdata=object$model,block.size=chunk.size,cluster=cluster))
} else { ## store discretization specific information to help with discrete prediction
  #object$dinfo <- list(gp=gp, v = G$v, ts = G$ts, dt = G$dt, qc = G$qc, drop = G$drop, pmf.names=pmf.names,lpip=lpip)
  #if (paratens) object$dinfo$para.discrete <- TRUE
  object$dinfo <- G$dinfo
}
rm(G);if (gc.level>0) gc()

if (is.null(object$fitted.values)) object$fitted.values <- family$linkinv(object$linear.predictors)

object$residuals <- if (is.null(family$residuals)) sqrt(family$dev.resids(object$y,object$fitted.values,object$prior.weights)) *
  sign(object$y-object$fitted.values) else residuals(object)
if (rho!=0) object$std.rsd <- AR.resid(object$residuals,rho,object$model$"(AR.start)")

if (!efam || is.null(object$deviance)) object$deviance <- sum(object$residuals^2)
## 'dev' is used in family$aic to estimate scale. That's standard and fine for Gaussian data, but
## can lead to badly biased estimates for e.g. low count data with the Tweedie (see Fletcher Biometrika paper)
## So set dev to give object $sig2 estimate when divided by sum(prior.weights)...
dev <- if (family$family!="gaussian"&&!is.null(object$sig2)) object$sig2*sum(object$prior.weights) else object$deviance ## used to give scale in family$aic
if (rho!=0&&family$family=="gaussian") dev <- sum(object$std.rsd^2)
object$aic <- if (efam) family$aic(object$y,object$fitted.values,family$getTheta(),object$prior.weights,dev) else
  family$aic(object$y,1,object$fitted.values,object$prior.weights,dev)
object$aic <- object$aic -
  2 * (length(object$y) - sum(sum(object$model[["(AR.start)"]])))*log(1/sqrt(1-rho^2)) + ## correction for AR
  2*sum(object$edf)
if (!is.null(object$edf2)&&sum(object$edf2)>sum(object$edf1)) object$edf2 <- object$edf1
if (is.null(object$null.deviance)) object$null.deviance <- sum(family$dev.resids(object$y,weighted.mean(object$y,object$prior.weights),object$prior.weights))
if (!is.null(object$full.sp)) {
  if (length(object$full.sp)==length(object$sp)&&
      all.equal(object$sp,object$full.sp)==TRUE) object$full.sp <- NULL
}
environment(object$formula) <- environment(object$pred.formula) <-
  environment(object$terms) <- environment(object$pterms) <-
  environment(attr(object$model,"terms"))  <- .GlobalEnv
if (control$trace) {
  t5 <- proc.time()
  t5 <- rbind(t1-t0,t2-t1,t3-t2,t4-t3,t5-t4)[,1:3]
  row.names(t5) <- c("initial","gam.setup","pre-fit","fit","finalise")
  print(t5)
}
names(object$gcv.ubre) <- method
object
} ## end of bam
