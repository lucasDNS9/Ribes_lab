#Hypergeometric test#

# Define parameters
q <- 2 #size of overlap 
m <- 174 #nb of gene in the gene set
n <- 11956 - m #nb of gene in total - nb of gene in the gene set
k <- 30 #nb of gene in the cluster

# Calculate p-value using hypergeometric distribution
p_value <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
#lower.tail True means p(X<k), if False: p(X>k)

# Print p-value
print(p_value)
