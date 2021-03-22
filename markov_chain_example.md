Markov Chain Example
================
Marco Smolla
3/22/2021

## What is a Markov Chain?

Markov Chains describe the transition probabilities of a system from one
state to another. For example, if a system has the state ‘activated’,
‘stand-by’, and ‘deactivated’, what is the probability to go from
‘stand-by’ to ‘activated’, and so on. For more information, make sure
to take a look at the corresponding [Wikipedia
entry](https://en.wikipedia.org/wiki/Markov_chain) and these [lecture
notes](https://www.stat.auckland.ac.nz/~fewster/325/notes/ch8.pdf).

## Which packages do we need?

There are a bunch of good `R` packages out there. Here, I am relying on
`markovchain` and for plotting I will use `igraph`.

``` r
# For calculating markov chains 
library(markovchain) 
```

    ## Package:  markovchain
    ## Version:  0.8.5-4
    ## Date:     2021-01-07
    ## BugReport: https://github.com/spedygiorgio/markovchain/issues

``` r
# For plotting as a network
library(igraph) 
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

## Example data

Here is our example data:

``` r
# Example transitions for a single individual (H=helper, R=reciprocator, B=beneficiary)
x <- sample(c("H", "R", "B"), 1000, replace = T, prob=c(.30, .60, .10)) # here I am creating a random vector of states (which in your case would be the sequence of states of your individual buyers)

# Count: 
table(x)
```

    ## x
    ##   B   H   R 
    ##  84 319 597

``` r
# Proportion
table(x)/length(x)
```

    ## x
    ##     B     H     R 
    ## 0.084 0.319 0.597

## Calculating a transition matrix - the packaged way

If you prefer to use a ready-made function, use the `markovchainFit()`
function from the `markovchain` package:

``` r
# using markovchainFit from the markovchain package
mcFit <- markovchainFit(data=x)
# look at results
show(mcFit)
```

    ## $estimate
    ## MLE Fit 
    ##  A  3 - dimensional discrete Markov Chain defined by the following states: 
    ##  B, H, R 
    ##  The transition matrix  (by rows)  is defined as follows: 
    ##            B         H         R
    ## B 0.04761905 0.3333333 0.6190476
    ## H 0.09717868 0.2821317 0.6206897
    ## R 0.08221477 0.3372483 0.5805369
    ## 
    ## 
    ## $standardError
    ##            B          H          R
    ## B 0.02380952 0.06299408 0.08584646
    ## H 0.01745381 0.02973929 0.04411049
    ## R 0.01174497 0.02378766 0.03120986
    ## 
    ## $confidenceLevel
    ## [1] 0.95
    ## 
    ## $lowerEndpointMatrix
    ##              B         H         R
    ## B 0.0009532288 0.2098672 0.4507916
    ## H 0.0629698436 0.2238437 0.5342347
    ## R 0.0591950491 0.2906254 0.5193667
    ## 
    ## $upperEndpointMatrix
    ##            B         H         R
    ## B 0.09428487 0.4567995 0.7873036
    ## H 0.13138752 0.3404196 0.7071447
    ## R 0.10523448 0.3838713 0.6417071
    ## 
    ## $logLikelihood
    ## [1] -877.5118

# Calculating a transition matrix - DIY

Alternatively, we can build our own `trans.matrix()` function:

``` r
trans.matrix <- function(X, prob=T)
{
 # turn sequence into a factor 
 if(!is.factor(x)) x <- as.factor(x)
 # define individual levels of the factor
 ids <- levels(x)
 # create an empty ids x ids matrix
 m <- matrix(0, ncol=length(ids), nrow=length(ids))
 # turn sequence into levle numbers (1, 2, 3, ...)
 state <- as.numeric(x)
 # create a two column matrix where the first column are t and the second column t+1 entries 
 tmp <- cbind(state[-length(state)], state[-1])
 # add up all individual transitions (i.e. at t the state is 1 and at t+1 the state is 2, then add 1 to the number in the matrix M_12)
 for(i in 1:nrow(tmp)){
  m[tmp[i,1], tmp[i,2]] <- m[tmp[i,1], tmp[i,2]] + 1
 }
 # turn counts into probabilities
 if(prob) m <- m / rowSums(m)
 # turn matrix into data.frame with column names and row names 
 df <- data.frame(m, row.names = ids)
 colnames(df) <- ids
 # return transitio matrix
 return(df)
}
```

Let us test this function without initial sequence `x`:

``` r
# calculate transition matrix
mcFit2 <- trans.matrix(x)
mcFit2
```

    ##            B         H         R
    ## B 0.04761905 0.3333333 0.6190476
    ## H 0.09717868 0.2821317 0.6206897
    ## R 0.08221477 0.3372483 0.5805369

Compare `mcFit2` with the estimate of `mcFit` and you’ll see that this
produces the same result. Of course, when using the markovchain package
you receive much more information than the transition matrix.

## Plot transition matrix as a graph

Finally, let us plot the result as an `igraph` network:

``` r
# turn transition matrix into a network
net <- graph_from_adjacency_matrix(adjmatrix = as.matrix(mcFit2), weighted = T)
# plot transiton matrix as network
plot(net, 
     edge.width=E(net)$weight*5, # use transition probabilities to vary edge widths
     edge.curved=T, # make edges curved (instead of straight)
     # vertex.size=50, # this is for fixed vertex size
     vertex.size=colSums(mcFit2)*50, # use in-edge weights (i.e. column sums of the transition matrix) to set the size of a vertex (larger vertices have larger weights of incoming edges)
     edge.label=round(E(net)$weight, 2)) # add transition probabilities as labels to all edges (and round values to make it easier to read)
```

![](markov_chain_example_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
