library(GeneCoexp)

# single
res <- nrcor(iris[,5],t(iris[,1:4]))

# multi
obj <- multinrcor(t(iris[,1:4]))

# create links
links <- create_links(obj)

# create network
net <- create_network(obj)
