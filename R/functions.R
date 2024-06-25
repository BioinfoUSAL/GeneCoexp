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

insidenrcor <- function(x, y, r.value, cor_function, sample.size, iter){
    matresults <- matrix(0,ncol=2,nrow=nrow(y))
    rownames(matresults) <- rownames(y)
    colnames(matresults) <- c("N","R")
    sd0<-apply(y, 1, sd)<0.1

    for(k in seq_len(iter)){
        msamp <- y[,sample(ncol(y), sample.size, TRUE),drop=FALSE]
        vsamp <- x[colnames(msamp)]
        aaa <- cor_function(msamp, vsamp)
        aaa[is.na(aaa)]<-0
        aaa[sd0,]<-0
        matresults[,"N"] <- matresults[,"N"] + (abs(aaa)>r.value)
        matresults[,"R"] <- matresults[,"R"] + (aaa)
    }
    #matresults[order(matresults[,1], decreasing = TRUE),, drop=FALSE]
    return(matresults)
}

nrcor <- function(x, y, r.value=NA, cor.method = c("pearson", "spearman"),
        sample.size=10, iter=10000){
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
        sample.size,iter)

    nmean <- results[,'N']
    pvalue <- pnorm(nmean, mean(nmean), sd(nmean), lower.tail=FALSE)
    results <- cbind(results,pvalue)
    padjust <- p.adjust(pvalue, method="fdr")
    results <- cbind(results,padjust)

    return(results)
}

nrcorcall <- function(ind, data, sample.size, cor.method, iter){
    if(ind<nrow(data)){
        vector <- data[ind,]
        matriz <- data[(ind+1):nrow(data),]
        if(is.matrix(matriz) == FALSE){
            matriz <- t(matriz)
            rownames(matriz) <- rownames(data)[(ind+1):nrow(data)]
        }
        r.value <- r.val(sample.size)
        insidenrcor(vector, matriz, r.value, get_cor_function(cor.method),
            sample.size, iter)
    }
}

multinrcor <- function(data, cor.method = c("pearson", "spearman"),
        sample.size = 10, iter = 10000, threads = NULL){
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
    clusterExport(cl, c("nrcor", "r.val", "cor.pearson", "cor.spearman"),
        envir=environment())

    results <- parLapply(cl, seq_len(nrow(data)-1), nrcorcall, data,
        sample.size, cor.method, iter)
    names(results) <- rownames(data)[-nrow(data)]
    stopCluster(cl)

    # Initialize a matrix to store the N and R columns
    N_columns <- R_columns <- matrix(0, nrow = nrow(data),
        ncol = nrow(data), dimnames = list(rownames(data),
        rownames(data)))
    # Fill the matrix with the N and R columns,
    for (r in seq_along(results)) {
        df <- results[[r]]
        N_columns[seq_len(nrow(df))+r, r] <- df[, 1]
        R_columns[seq_len(nrow(df))+r, r] <- df[, 2]
    }
    N_columns[upper.tri(N_columns)] <- t(N_columns)[upper.tri(N_columns)]
    R_columns[upper.tri(R_columns)] <- t(R_columns)[upper.tri(R_columns)]

    pv <- multinrcor_p(N_columns)

    structure(list(N = N_columns, R = R_columns, pvalue = pv[[1]],
        padjust = pv[[2]], iter = iter), class = "multinrcor")
}

multinrcor_p <- function(N_columns){
    nvalues <- as.vector(N_columns)
    means <- as.vector(rep(rowMeans(N_columns), each=nrow(N_columns)))
    mads <- as.vector(rep(apply(N_columns, 1, mad), each=nrow(N_columns)))

    pvalue <- pnorm(nvalues, means, mads, lower.tail=FALSE)
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

    return(list(pvalue_columns,padjust_columns))
}

create_links <- function(x, cutoff = NULL){
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
        linksfilter <- links[links$padjust<=cutoff,]
    }else{
        linksfilter <- links
    }

    return(linksfilter)
}

create_network <- function(x, cutoff = NULL, nodes = NULL, name = NULL,
        label = NULL){
    if(!inherits(x,'multinrcor')){
        stop("x: must be a multinrcor object")
    }
    if(!requireNamespace("rD3plot", quietly=TRUE)){
        stop("Install 'rD3plot' to create networks.")
    }

    links <- create_links(x,cutoff)
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
    return(net)
}
