r.val <- function(n, alpha=.01){
    df <- n-2
    critical.t <- qt(alpha/2, df, lower.tail = FALSE)
    critical.r <- sqrt((critical.t^2)/((critical.t^2)+df))
    return(critical.r)
}

cor.pearson <- function(msamp, vsamp){
    mmean <- rowMeans(msamp)
    vmean <- mean(vsamp)
    n<-length(vsamp)
    covar <- ((msamp - mmean)%*%(vsamp - vmean))/(n-1)
    sd1 <- matrixStats::rowSds(msamp)
    sd2 <- sqrt(sum((vsamp-vmean)^2)/(n-1))
    covar/(sd1 * sd2)
}

cor.spearman <- function(msamp, vsamp){
    mrank <- t(matrixStats::rowRanks(msamp, ties.method="average"))
    vrank <- rank(vsamp)
    mrankmean <- colMeans(mrank)
    vrankmean <- mean(vrank)
    num <- t(mrank - mrankmean) %*% (vrank - vrankmean)
    denom <- sqrt(colSums((mrank - mrankmean)^2) * sum((vrank - vrankmean)^2))
    num/denom
}

get_cor_function <- function(cor.method){
    cor.method <- cor.method[1]
    cor_function <- NULL
    if(cor.method == "pearson"){
        cor_function <- cor.pearson
    }else if(cor.method == "spearman"){
        cor_function <- cor.spearman
    }else{
        stop("cor.method: invalid correlation method")
    }
    return(cor_function)
}

insidenrcor <- function(x, y, r.value, cor_function, sample.size, iter,
        samples_results){
    resultsN <- vector(mode="numeric", length=nrow(y))
    resultsR <- vector(mode="numeric", length=nrow(y))
    names(resultsN) <- names(resultsR) <- rownames(y)

    if(samples_results){
        samplesmat <- matrix(0, nrow = ncol(y), ncol = ncol(y),
            dimnames = list(colnames(y), colnames(y)))
        samplesnum <- vector(mode="numeric", length=ncol(y))
        names(samplesnum) <- colnames(y)
    }else{
        samplesmat <- NULL
        samplesnum <- NULL
    }

    sd0<-apply(y, 1, sd)<0.1

    for(k in seq_len(iter)){
        msamp <- y[,sample(ncol(y), sample.size, TRUE),drop=FALSE]
        vsamp <- x[colnames(msamp)]
        cormat <- cor_function(msamp, vsamp)
        cormat[is.na(cormat)]<-0
        cormat[sd0,]<-0
        abscormatgtrvalue <- abs(cormat)>r.value

        resultsN <- resultsN + abscormatgtrvalue
        resultsR <- resultsR + (cormat*abscormatgtrvalue)
        if(samples_results){
            samplesmat[colnames(msamp),colnames(msamp)] <-
        samplesmat[colnames(msamp),colnames(msamp)] + sum(abscormatgtrvalue)
            samplesnum[colnames(msamp)] <- samplesnum[colnames(msamp)]+1
        }
    }
    return(list(N = resultsN, R = resultsR,
        samplesmat = samplesmat, samplesnum = samplesnum))
}

create_structure <- function(thisclass, N, R, pvalue, padjust,
        pvaluev, padjustv, samplesmat, samplesnum, iter){
    l <- list(N = N, R = R, pvalue = pvalue, padjust = padjust)
    if(!is.null(pvaluev)){
        l[['pvaluev']] <- pvaluev
    }
    if(!is.null(padjustv)){
        l[['padjustv']] <- padjustv
    }
    if(!is.null(samplesmat)){
        l[['samplesmat']] <- samplesmat
    }
    if(!is.null(samplesnum)){
        l[['samplesnum']] <- samplesnum
    }
    l[['iter']] <- iter
    structure(l, class = thisclass)
}

nrcor <- function(x, y, r.value=NA, cor.method = c("pearson", "spearman"),
        mads = FALSE, sample.size=10, iter=10000, samples_results=TRUE){
    if(!is.numeric(x)){
        x <- as.numeric(x)
    }

    if(!is.matrix(y)){
        y <- data.matrix(y)
    }

    if(is.null(colnames(y))){
        colnames(y) <- paste0("sample_",seq_len(ncol(y)))
    }

    if(length(x)!=ncol(y)){
        stop("'x' length must match with 'y' number of columns")
    }

    if(is.null(names(x))){
        names(x) <- colnames(y)
    }else if(length(setdiff(colnames(y),names(x)))){
        stop("'x' names must match with 'y' column names")
    }

    if(is.na(r.value)){
        r.value <- r.val(sample.size)
    }

    results <- insidenrcor(x,y,r.value,get_cor_function(cor.method),
        sample.size,iter,samples_results)

    nvalues <- results[['N']]
    means <- mean(nvalues)
    if(mads){
        pnormsd <- mad(nvalues)
    }else{
        pnormsd <- sd(nvalues)
    }
    pvalue <- pnorm(nvalues, means, pnormsd, lower.tail=FALSE)
    padjust <- p.adjust(pvalue, method="fdr")

    return(create_structure("nrcor", nvalues, results[['R']],
        pvalue, padjust, NULL, NULL,
        results[['samplesmat']], results[['samplesnum']], iter))
}

nrcorcall <- function(ind, data, sample.size, cor.method, iter,
        samples_results){
    if(ind<nrow(data)){
        vector <- data[ind,]
        matriz <- data[(ind+1):nrow(data),]
        if(is.matrix(matriz) == FALSE){
            matriz <- t(matriz)
            rownames(matriz) <- rownames(data)[(ind+1):nrow(data)]
        }
        r.value <- r.val(sample.size)
        insidenrcor(vector, matriz, r.value, get_cor_function(cor.method),
            sample.size, iter, samples_results)
    }
}

multinrcor <- function(data, cor.method = c("pearson", "spearman"),
        mads = FALSE, pvalue_var = FALSE, sample.size = 10, iter = 10000,
        samples_results = TRUE, threads = NULL){
    if(!is.matrix(data)){
        data <- data.matrix(data)
    }
    if(is.null(colnames(data))){
        colnames(data) <- paste0("sample_",seq_len(ncol(data)))
    }
    if(!is.numeric(threads)){
        threads <- detectCores()
    }
    cl <- makeCluster(threads)

    results <- parLapply(cl, seq_len(nrow(data)-1), nrcorcall, data,
        sample.size, cor.method, iter, samples_results)
    names(results) <- rownames(data)[-nrow(data)]
    stopCluster(cl)

    # Initialize a matrix to store the N and R columns
    N_columns <- R_columns <- matrix(0, nrow = nrow(data),
        ncol = nrow(data), dimnames = list(rownames(data),
        rownames(data)))
    if(samples_results){
        samplesmat <- matrix(0, nrow = ncol(data), ncol = ncol(data),
            dimnames = list(colnames(data), colnames(data)))
        samplesnum <- vector(mode="numeric", length=ncol(data))
    }else{
        samplesmat <- samplesnum <- NULL
    }
    # Fill the matrix with the N and R columns,
    for(r in seq_along(results)) {
        resN <- results[[r]][['N']]
        resR <- results[[r]][['R']]
        N_columns[seq_along(resN)+r, r] <- resN
        R_columns[seq_along(resR)+r, r] <- resR

        if(samples_results){
            samplesmat <- samplesmat + results[[r]][['samplesmat']]
            samplesnum <- samplesnum + results[[r]][['samplesnum']]
        }
    }
    N_columns[upper.tri(N_columns)] <- t(N_columns)[upper.tri(N_columns)]
    R_columns[upper.tri(R_columns)] <- t(R_columns)[upper.tri(R_columns)]

    nvalues <- as.vector(N_columns)
    means <- mean(nvalues)
    pnormsdfn <- sd
    if(mads){
        pnormsdfn <- mad
    }
    pnormsd <- pnormsdfn(nvalues)
    pvalue <- round(pnorm(nvalues, means, pnormsd, lower.tail=FALSE),7)
    pvalue_columns <- matrix(pvalue, nrow = nrow(N_columns))

    for(i in seq_len(ncol(pvalue_columns))){
        for(j in seq_len(nrow(pvalue_columns))){
            pvalue_columns[i,j] <- pvalue_columns[j,i] <-
                min(pvalue_columns[i,j], pvalue_columns[j,i])
        }
    }

    padjust <- p.adjust(pvalue, method="fdr")
    padjust_columns <- matrix(padjust, nrow = nrow(N_columns))

    for(i in seq_len(ncol(padjust_columns))){
        for(j in seq_len(nrow(padjust_columns))){
            padjust_columns[i,j] <- padjust_columns[j,i] <-
                min(padjust_columns[i,j], padjust_columns[j,i])
        }
    }

    diag(pvalue_columns) <- diag(padjust_columns) <- 1
    rownames(pvalue_columns) <- rownames(padjust_columns) <- rownames(N_columns)
    colnames(pvalue_columns) <- colnames(padjust_columns) <- rownames(N_columns)

    pvalue_columnsv <- NULL
    padjust_columnsv <- NULL
    if(pvalue_var){
        meansv <- as.vector(rep(rowMeans(N_columns), times=nrow(N_columns)))
        pnormsdfnv <- sd
        if(mads){
            pnormsdfnv <- mad
        }
        pnormsdv <- as.vector(rep(apply(N_columns, 1, pnormsdfnv),
            times=nrow(N_columns)))

        pvaluev <- round(pnorm(nvalues, meansv, pnormsdv, lower.tail=FALSE),7)
        pvalue_columnsv <- matrix(pvaluev, nrow = nrow(N_columns))
    
        for(i in seq_len(ncol(pvalue_columnsv))){
            for(j in seq_len(nrow(pvalue_columnsv))){
                pvalue_columnsv[i,j] <- pvalue_columnsv[j,i] <-
                    min(pvalue_columnsv[i,j], pvalue_columnsv[j,i])
            }
        }
    
        padjustv <- p.adjust(pvaluev, method="fdr")
        padjust_columnsv <- matrix(padjustv, nrow = nrow(N_columns))
    
        for(i in seq_len(ncol(padjust_columnsv))){
            for(j in seq_len(nrow(padjust_columnsv))){
                padjust_columnsv[i,j] <- padjust_columnsv[j,i] <-
                    min(padjust_columnsv[i,j], padjust_columnsv[j,i])
            }
        }
    
        diag(pvalue_columnsv) <- diag(padjust_columnsv) <- 1
        rownames(pvalue_columnsv) <- rownames(padjust_columnsv) <-
            rownames(N_columns)
        colnames(pvalue_columnsv) <- colnames(padjust_columnsv) <-
            rownames(N_columns)
    }

    return(create_structure("multinrcor", N_columns, R_columns,
        pvalue_columns, padjust_columns, 
        pvalue_columnsv, padjust_columnsv,
        samplesmat, samplesnum, iter))
}

create_links <- function(x, cutoff = NULL, cutoffvar = "padjust"){
    if(!inherits(x,'multinrcor')){
        stop("x: must be a multinrcor object")
    }

    N <- x$N
    R <- x$R
    pvalue <- x$pvalue
    padjust <- x$padjust
    source <- rep(rownames(N),nrow(N))
    target <- rep(colnames(N),rep(ncol(N),nrow(N)))
    links <- data.frame(source=source, target=target,
        N=as.vector(as.matrix(N)), R=as.vector(as.matrix(R)),
        pvalue=as.vector(as.matrix(pvalue)),
        padjust = as.vector(as.matrix(padjust)))
    links <- links[links$N>0,]

    if(is.numeric(cutoff)){
        if(!(cutoffvar %in% c("N","R","pvalue","padjust"))){
            warning("cutoffvar: must be 'N', 'R', 'pvalue' or 'padjust'")
            linksfilter <- links
        }else{
            linksfilter <- links[links[[cutoffvar]]<=cutoff,]
        }
    }else{
        linksfilter <- links
    }

    return(linksfilter)
}

create_network <- function(x, cutoff = NULL, cutoffvar = "padjust",
        nodes = NULL, name = NULL, label = NULL){
    if(!inherits(x,'multinrcor')){
        stop("x: must be a multinrcor object")
    }
    if(!requireNamespace("rD3plot", quietly=TRUE)){
        stop("Install 'rD3plot' to create networks.")
    }

    links <- create_links(x,cutoff,cutoffvar)
    if(nrow(links)>10000){
        warning(
"Too much links will cause performance issues in web browser.
You can use the cutoff."
        )
    }

    links[,'R/iter'] <- links[,'R'] / x$iter
    links[,'-log10pvalue'] <- -log10(links[,'pvalue'])

    if(!is.null(nodes)){
        if(is.data.frame(nodes)){
            if(is.null(names(nodes))){
                names(nodes) <- paste0("attr",seq_len(ncol(nodes)))
                if(!is.null(name) || !is.null(label)){
                    warning(
"You can't specify name or label if there are no nodes' colnames"
                    )
                    name <- NULL
                    label <- NULL
                }
            }
            if(is.null(name) || !(name %in% names(nodes))){
                name <- names(nodes)[1]
            }
            subnodes <- data.frame(attr1 = union(links$source,links$target))
            names(subnodes) <- name
            nodes <- merge(subnodes,nodes,by=name,all.x=TRUE)
        }else{
            warning("nodes: must be a data frame")
        }
    }

    net <- rD3plot::network_rd3(links=links, nodes=nodes, name=name,
        label=label, lcolor="R/iter", lwidth="-log10pvalue", linkBipolar=TRUE)
    plot(net)
    return(net)
}

plot.multinrcor <- plot.nrcor <- function(x, cutoff = 0.05,
        cutoffvar = "padjust", ...){
    colors <- rep("black",length(x$N))
    if(cutoffvar=="R"){
        colors[abs(x$R/x$N)<=cutoff] <- "red"
    }else{
        colors[x[[cutoffvar]]<=cutoff] <- "red"
    }
    plot(x$N, x$R/x$N, pch=16, col=colors)
}
