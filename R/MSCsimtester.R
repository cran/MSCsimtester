#' Validity tests of simulators of the multspecies coalescent model in phylogenomics.
#'
#' The package performs comparisons of certain summary statistics for simulated gene tree samples to
#' theoretical predictions under the multispecies coalescent model.
#' The primary functions are \code{rootedTriple} for comparison of frequencies of topological rooted triples
#' on gene trees, and \code{pairwiseDist} and \code{ADtest} for comparison of the distributions of pairwise
#' distances between taxa on gene trees.
#'
#' Required input is a collection of gene trees, stored as a \code{multiPhylo} object by the \code{ape} package,
#' and a rooted species tree, as a \code{Phylo} object,
#' with edge lengths in generations, together with constant population sizes for each edge.
#'
#' \code{MSCsimtester} builds on the packages \code{ape} and \code{kSamples}.
#'
#' For further examples of use and citation purposes, see \insertCite{testing19}{MSCsimtester}.
#'
#' @references
#' \insertRef{testing19}{MSCsimtester}
#'
#' @docType package
#' @name MSCsimtester
#'
#' @import stats ape kSamples
#' @importFrom graphics hist lines plot
#' @importFrom Rdpack reprompt
NULL


#' Plot species tree, with edge numbers on edges.
#'
#' Under the MSC, each edge in the species tree must be assigned a population
#' size.  This function displays the species tree with the edges
#' numbered, to aid the user in entering constant population sizes as an
#' appropriately ordered list.
#'
#' @param stree An object of class \code{phylo} containing a rooted metric
#' species tree.
#'
#' @examples stree=read.tree(text="(((a:10000,b:10000):10000,c:20000):10000,d:30000);")
#' plotEdgeOrder(stree)
#' pops=c(30000,20000,1,1,1,1,10000)
#' plotPops(stree,pops)
#'
#' @seealso \code{\link{pairwiseDist}}, \code{\link{rootedTriple}}, \code{\link{plotPops}}
#'
#' @return NONE
#'
#' @export
plotEdgeOrder = function(stree)
{
  if (is.rooted(stree) == FALSE) {
    stop("Tree should be rooted.")
  }
  stree$root.edge <- mean(stree$edge.length)
  plot.phylo(
    stree,
    label.offset = 1,
    cex = 1.3,
    edge.width = 3,
    main = "Edge Order",
    root.edge = TRUE
  )
  numEdges = dim(stree$edge)[1]
  edgeNumbers = c(1:numEdges)
  edgelabels(edgeNumbers,
             srt = 0,
             bg = "white",
             cex = 1.3)
  nodelabels(
    numEdges + 1,
    length(stree$tip.label) + 1,
    srt = 0,
    bg = "white",
    cex = 1.3,
    adj = 2
  )
}

#' Plot species tree, with population sizes on edges.
#'
#' @param stree An object of class \code{phylo} containing a rooted metric species tree.
#' @param populations A vector containing constant population sizes, one entry for each
#' edge/population in the species tree, with last entry for
#' the population ancestral to the root.
#'
#' @return NONE
#'
#' @examples stree=read.tree(text="(((a:10000,b:10000):10000,c:20000):10000,d:30000);")
#' plotEdgeOrder(stree)
#' pops=c(30000,20000,1,1,1,1,10000)
#' plotPops(stree,pops)
#'
#' @seealso \code{\link{pairwiseDist}}, \code{\link{rootedTriple}}, \code{\link{plotEdgeOrder}}
#'
#' @export
plotPops = function(stree, populations)
{
  if (is.rooted(stree) == FALSE) {
    stop("Tree should be rooted.")
  }
  len = length(populations)
  if (dim(stree$edge)[1] != len - 1) {
    stop("Number of population sizes is incorrect.")
  }
  rootPop = populations[len]
  populations = populations[1:len - 1]
  stree$root.edge <- mean(stree$edge.length)
  plot.phylo(
    stree,
    label.offset = 1,
    cex = 1.3,
    edge.width = 1.3,
    main = "Population Sizes",
    root.edge = TRUE
  )
  edgelabels(
    populations,
    srt = 0,
    bg = "white",
    cex = 1.5,
    col = "deepskyblue2"
  )
  nodelabels(
    rootPop,
    length(stree$tip.label) + 1,
    srt = 0,
    bg = "white",
    cex = 1.5,
    col = "deepskyblue2",
    adj = 1
  )
}


#' Compute and plot sample and theoretical pairwise distance densities.
#'
#' Computes theoretical pairwise distance densities under the MSC on a species tree and empirical pairwise distances from
#' gene trees in a sample. A histogram of empirical values is plotted over the theoretical pdf.
#'
#' @param stree An object of class \code{phylo} containing a rooted metric species tree. Edge lengths are in
#' number of generations.
#' @param popSizes A vector containing constant population sizes, one entry for each
#' edge/population in the species tree, for a haploid population.  Sizes should be doubled for diploids.
#' If \code{stree} has k edges, then \code{popSizes} must have k+1 elements, with final entry
#' the size of the population ancestral to the root.
#' @param gtSample An object of class \code{multiPhylo} holding a sample of gene trees from a simulation.
#' Taxon labels on gene trees must be identical to those on \code{stree}.
#' @param taxon1 A string specifying one taxon on \code{stree}.
#' @param taxon2 A string specifying a second taxon on \code{stree}, distinct
#' from \code{taxon1}.
#' @param numSteps A positive integer number of values to be computed for
#' graphing the theoretical pairwise distance density.  Default is \code{numSteps = 1000}.
#' A larger value produces a smoother plot.
#' @param tailProb A cutoff value, between 0 and 1, for the theoretical density, with a default of 0.01.
#' The theoretical pairwise distance is plotted from (0, xMax), where  xMax
#' is the larger of the maximum pairwise distance in the gene tree sample and the value cutting off
#' a tail of area \code{tailProb} under the pdf. A message returns the proportion of sample distances in this tail.
#'
#'
#' @return A list of items needed for Anderson-Darling test(s), for use by \code{ADtest},
#' returned invisibly. See function code for more details.
#'
#' @examples stree=read.tree(text="((((a:10000,b:10000):10000,c:20000):10000,d:30000):10000,e:40000);")
#' pops=c(15000,25000,10000,1,1,1,1,1,12000)
#' gts=read.tree(file=system.file("extdata","genetreeSample",package="MSCsimtester"))
#' pairwiseDist(stree,pops,gts,"a","b")
#'
#' @seealso \code{\link{plotEdgeOrder}}, \code{\link{plotPops}}, \code{\link{ADtest}}
#'
#' @export
#'
pairwiseDist = function(stree,
                        popSizes,
                        gtSample,
                        taxon1,
                        taxon2,
                        numSteps = 1000,
                        tailProb = .01)
{
  # some initial (limited) checks

  if (identical(taxon1, taxon2)) {
    warning('Taxa must be different. Exiting.')
    return(invisible())
  }

  # Check that the tip labels in species tree and the first gene tree agree.
  if (setequal(stree$tip.label, gtSample[[1]]$tip.label) == FALSE)  {
    warning('Tip labels on species tree and gene tree 1 do not match.  Exiting.')
    return(invisible())
  }

  # Check that a vector of numeric constant pop sizes was entered
  if (class(popSizes) != "numeric") {
    warning("Population sizes must be numeric. Exiting.")
    return(invisible())
  }

  N = length(gtSample)  # number of gene trees in sample

  # check timing for large gt datasets
  start_time <- Sys.time()

  # obtains distance in each gene tree betwen taxon1 and taxon2
  sampleDist = vector(mode = "numeric", length = N) # allocate space for sample distances
  for (z in 1:N)
  {
    tr = gtSample[[z]]  # get gene tree
    d = cophenetic(tr)  # compute distance table for given gene tree
    sampleDist[z] = d[taxon1, taxon2] # record the distance between the taxa of interest
  }
  current_time = Sys.time()
  elapsedTime = as.numeric(difftime(current_time, start_time, units = "mins"))
  if (elapsedTime > 4) {
    message(paste0(
      "  Time to read and process ",
      N,
      " gene tree samples was ",
      round(elapsedTime, 1),
      " minutes."
    ))
  }

  # create an empirical histogram for full sample for plotting later
  hg = hist(sampleDist, plot = FALSE, breaks = 100)

  # Create and plot theoretical pdf

  # invert the constant population sizes and save as functions
  # (overkill here but useful for nonconstant population sizes in future)
  invPopSizes = lapply(popSizes, invertPopSize)

  edgeLens = stree$edge.length

  # holds the edge numbers of the edges above taxon1 and taxon2
  ancestry = edgesabove(stree, taxon1, taxon2)
  ea1 = edgesabove(stree, taxon1, taxon1)
  ea2 = edgesabove(stree, taxon2, taxon2)

  # get edge numbers for the edges from taxon i to MRCA
  t1 = ea1[!ea1 %in% ancestry]
  t2 = ea2[!ea2 %in% ancestry]

  # compute the number of internal nodes above MRCA of taxon1, taxon2
  numInternalNodes = length(ancestry - 1)

  # sum distances from internal nodes to two tips taxon1 and taxon2
  # entries will correspond to g_ab + 2*m_i in paper
  distancesAtInternalNodes = vector(mode = "numeric", length = numInternalNodes)
  distancesAtInternalNodes[1] = sum(stree$edge.length[t1], stree$edge.length[t2])
  if (numInternalNodes > 1) {
    for (i in 2:numInternalNodes) {
      distancesAtInternalNodes[i] = distancesAtInternalNodes[i - 1] + 2 * edgeLens[ancestry[i -
                                                                                              1]]
    }
  }

  # vector to hold the probabilities of no coalescence occuring before reaching
  # internal nodes. These are the products of what are called eta in paper.
  nocoalProbAtInternalNodes = vector(mode = "numeric", length = numInternalNodes)
  nocoalProbAtInternalNodes[1] = 1
  if (numInternalNodes > 1) {
    for (i in 2:numInternalNodes) {
      nocoalProbAtInternalNodes[i] = nocoalProbAtInternalNodes[i - 1] *
        exp(-edgeLens[ancestry[i - 1]] * invPopSizes[[ancestry[i - 1]]](0))
    }
  }

  # get the x values for plotting the theoretical distribution

  # get the minimal distance with non-zero support
  # do not get maximum distance since later this code will be adapted for non-ultrametric trees
  d = cophenetic(stree)
  xMin = d[taxon1, taxon2]

  # get xMax based on theoretical distribution and sample size N
  indRootPop = ancestry[numInternalNodes]

  # get x value that cuts off a tail of tailProb
  xMax = round(
    distancesAtInternalNodes[numInternalNodes] -
      2 * log(tailProb / nocoalProbAtInternalNodes[numInternalNodes]) / (invPopSizes[[indRootPop]](0))
  )
  message(
    paste0(
      "  Proportion of sample distances in the ", tailProb,"-tail of theoretical distribution is ",round(sum(sampleDist > xMax) / N, 3),"."
    )
  )

  xMax = max(xMax, ceiling(max(sampleDist)))

  xStepLength = (xMax - xMin) / (numSteps - 1)

  xVals = seq(xMin, xMax, xStepLength)      # the domain of the pdf

  # allocate space for the theoretical density
  theorDensity <- vector(mode = "numeric", length = numSteps)

  # compute the theoretical density

  inds = findInterval(xVals, distancesAtInternalNodes)

  numIntervals = inds[numSteps]

  for (i in 1:numIntervals) {
    xx = xVals[inds == i]
    theorDensity[inds == i] = nocoalProbAtInternalNodes[i] *
      (exp(-(xx - distancesAtInternalNodes[i]) * invPopSizes[[ancestry[i]]](0) /
             2)) *
      (invPopSizes[[ancestry[i]]](0)) / 2

  }

  # Now plot both theoretical density and sample distances histogram
  title = paste0("Density of d(", taxon1, ",", taxon2, ") on gene trees")

  # plot theoretical density first
  plot(
    c(0, xVals[1], xVals),
    c(0, 0, theorDensity),
    col = "blue",
    lwd = 2,
    type = 'l',
    xlab = "Generations",
    ylab = "density",
    main = title,
    ylim = c(0, max(theorDensity, hg$density))
  )

  # plot empirical histogram for full sample
  plot(hg,
       xlab = "Generations",
       freq = FALSE,
       add = TRUE)

  return(invisible(
    list(
      sampleDist,
      distancesAtInternalNodes,
      nocoalProbAtInternalNodes,
      xMin,
      xMax,
      max(theorDensity),
      taxon1,
      taxon2,
      ancestry,
      invPopSizes
    )
  ))
  # The returned list here contains information to be passed to the ADtest function if a formal
  # statistical test of the distribution fit is desired. These include
  # sampleDist: the vector of distances between taxon1 and taxon 2 from the sample
  # distancesAtInternalNodes: the sums of the 2 distances from taxon1 and taxon2 to ancestral internal nodes
  # nocoalProbAtInternalNodes: the probabilities of no coalescence between the two lineages below internal nodes
  # xMin: the smallest distance at which a coalescent event could occur
  # xMax: the largest distance which will be included in plot
  # max(theorDensity): The maximum of the theoretical density
  # ancestry: a vector of the edges ancestral to taxon1 and taxon2
}


#' Anderson-Darling test comparing sample and theoretical pairwise distance distributions.
#'
#' Takes as input theoretical pairwise distance densities under the MSC and
#' empirical pairwise distances from gene trees in a sample, as returned by
#' the function \code{pairwiseDist}.  Uses the package \code{kSamples} to perform
#' either one test on the entire dataset or multiple tests on subsamples.
#'
#' @details The Anderson-Darling test compares the empirical distance distribution for a supplied gene tree
#' sample to a sample drawn from the theoretical distribution. The output, passed from the \code{kSamples} package,
#' thus says that 2 samples are being compared, to test a null-hypothesis that they come from the same distribution.
#' See \code{kSamples} documentation for function \code{ad.test} for more details.
#'
#' Repeated runs of this function will give different results, since the sample from the theoretical distribution
#' will vary. Under the null hypothesis p-values for different runs should be approximately
#' uniformly distributed.
#'
#' Numerical issues may result in poor performance of Anderson-Darling tests when the sample size
#' is very large, so
#' an optional parameter \code{subsampleSize} can be set to create subsamples of smaller size.
#' If \code{subsampleSize} is a positive integer,
#' Anderson-Darling tests are performed on each subset, comparing them to
#' a random sample of the same size from the theoretical distribution. Good fit is indicated by an approximately uniform
#' distribution of the subsample p-values.
#'

#'
#' @param distanceDensities A list containing values needed for performing Anderson-Darling
#' test(s) on a gene tree sample and species tree, as output by \code{pairwiseDist}.
#' For details, see code for \code{pairwiseDist}.
#'
#' @param subsampleSize A positive integer to perform multiple tests on subsamples,
#' or \code{FALSE} (default) to perform one test on full sample.
#'
#' @return An object of type \code{ADtestOutput} including a sample \code{$Sample} from the theoretical distance distribution of
#'  the same size as the empirical one, and \code{$ADtest} which is of type \code{kSamples} and
#'  has all output from the Anderson-Darling test if only
#'  one test was performed, or the number of tests if tests were performed on subsamples.
#'
#' @examples \donttest{stree=read.tree(text="((((a:10000,b:10000):10000,c:20000):10000,d:30000):10000,e:40000);")
#' pops=c(15000,25000,10000,1,1,1,1,1,12000)
#' gts=read.tree(file=system.file("extdata","genetreeSample",package="MSCsimtester"))
#' distDen=pairwiseDist(stree,pops,gts,"a","b")
#' ADtest(distDen)
#' ADtest(distDen,1000) }
#'
#' @seealso \code{\link{pairwiseDist}}, \code{\link{kSamples-package}}
#'
#' @export
#'
ADtest = function(distanceDensities, subsampleSize = FALSE)
{
  # returned invisibly by pairwiseDist()
  # return(invisible(list(sampleDist,distancesAtInternalNodes,nocoalProbAtInternalNodes,
  #                       xMin,xMax,max(theorDensity),taxon1,taxon2,ancestry,invPopSizes)))

  empirSample = distanceDensities[[1]]
  distancesAtInternalNodes = distanceDensities[[2]]
  nocoalProbAtInternalNodes = distanceDensities[[3]]
  xMin = distanceDensities[[4]]
  xMax = distanceDensities[[5]]
  C = distanceDensities[[6]]
  taxon1 = distanceDensities[[7]]
  taxon2 = distanceDensities[[8]]
  ancestry = distanceDensities[[9]]
  invPopSizes = distanceDensities[[10]]

  # length of empircal density (= number of gene trees in sample)
  N = length(empirSample)

  numInternalNodes = length(distancesAtInternalNodes)
  distanceAtRoot = distancesAtInternalNodes[numInternalNodes]

  if (subsampleSize == FALSE) {
    subsampleSize = N
  }
  if ((subsampleSize <= 0) || (subsampleSize > N)) {
    warning('Test size for AD test must be positive and less than the number of gene trees.  Exiting.')
    return(invisible())
  }

  # perform Anderson-Darling tests
  # compute number of AD tests to perform
  if (subsampleSize == N) {
    numTests = 1
  }
  else {
    numTests = floor(N / subsampleSize)
  }

  # create sample from theoretical distribution
  theorSample <- vector(mode = "numeric", length = N)

  start_time = Sys.time()

  k <- 1 # sample counter
  while (k < N + 1)
  {
    z <- runif(n = 1, min = xMin, max = xMax)

    # compute pairwise density for z
    if (z >= distanceAtRoot) {
      ind = numInternalNodes
    }
    else {
      ind = which.min(z >= distancesAtInternalNodes) - 1
    }
    R = nocoalProbAtInternalNodes[ind] *
      (exp(-(z - distancesAtInternalNodes[ind]) * invPopSizes[[ancestry[ind]]](0) /
             2) *
         invPopSizes[[ancestry[ind]]](0) / 2)  / C

    if (R > runif(1))
    {
      theorSample[k] <- z
      k  <- k + 1
    }
  }

  current_time = Sys.time()
  elapsedTime = as.numeric(difftime(current_time, start_time, units = "mins"))
  if (elapsedTime > 4) {
    message(
      paste0(
        "  Total time to compute ",
        k - 1,
        " samples from theoretical density was ",
        round(elapsedTime, 1),
        " minutes."
      )
    )
  }

  pval = rep(0, numTests)  # variable to hold the p-values
  for (j in 1:numTests) {
    # perform numTests ad.tests and save p-value
    subind = (1:subsampleSize) + (j - 1) * subsampleSize
    s = ad.test(empirSample[subind], theorSample[subind])
    pval[j] = s$ad[1, 3]
  }
  if (numTests == 1) {
    obj = list()
    obj$Adtest = s
    obj$Sample = invisible(theorSample)
    class(obj) = "ADtestOutput"
  }
  else {
    hist(
      pval,
      breaks = numTests,
      main = "Anderson-Darling test results",
      xlab = "p-values",
      xlim = c(0, 1),
      freq = TRUE
    )
    obj = list()
    obj$Adtest = numTests
    obj$Sample = invisible(theorSample)
    class(obj) = "ADtestOutput"
  }
  return(obj)
}




#' Compare expected and sample frequencies of topological rooted triples.
#'
#' For a given species tree with population sizes, compares the expected frequencies
#' of rooted triples to empirical frequencies in a sample of gene trees, using Chi-squared tests with 2 d.f.
#' The exact and estimated internal branch length (in coalescent units) of the rooted triple
#' in the species tree are also computed for comparison.  A single test can be performed
#' on the entire gene tree sample, or multiple tests on subsamples.
#'
#' @details
#' When \code{subsampleSize} is \code{FALSE} the Chi-squared test is performed using all
#' gene trees in \code{gtSample}. Results are reported in tabular form in the console.
#'
#' @details
#' When \code{subsampleSize} is positive, the  \code{N} trees in \code{gtSample}
#' will be partitioned into \code{N/subsampleSize} subsamples, with a Chi-squared test
#' performed for each. Histograms are plotted for (1) the p-values for the Chi-squared tests on
#' subsamples, and (2) subsample estimates of the internal branch length for the rooted triple on the
#' species tree, with the true value marked.
#'
#' @details Three distinct taxon names must be supplied, all of which must occur on \code{stree}
#' and in each of the gene trees in the sample.
#'
#' @param stree An object of class \code{phylo} containing a rooted metric species tree.
#' Edge lengths are in number of generations.
#' @param popSizes An ordered list containing constant population sizes for each species tree edge, for a haploid organism.
#' Sizes should be doubled for diploids. If \code{stree} has k edges, then \code{popSizes} must have k+1 elements,
#' with the final entry for the population ancestral to the root.
#' @param gtSample An object of class \code{multiPhylo} holding a sample of gene trees from a simulation.
#' Taxon labels on gene trees must be identical to those on \code{stree}.
#' @param taxon1 A string specifying one taxon on \code{stree}.
#' @param taxon2 A string specifying a second taxon on \code{stree}, distinct from \code{taxon1}.
#' @param taxon3 A string specifying a third taxon on \code{stree}, distinct from \code{taxon1}, \code{taxon2}.
#' @param subsampleSize A positive integer or \code{FALSE}, giving size of subsamples of \code{gtSample} to analyze.
#'
#' @return If \code{subsampleSize} is \code{FALSE}, returns an object of type \code{rootedTripleOutput}
#' which contains a table \code{$TripletCounts} of empirical and expected rooted triple counts,
#' a p-value \code{$pv} from the Chi-squared test, and a column \code{$InternalEdge} of estimated and exact internal edge lengths.
#' If \code{subsampleSize} is \code{TRUE}, returns \code{NULL} but produces several plots.
#'
#' @seealso \code{\link{plotEdgeOrder}}, \code{\link{plotPops}}
#' @examples \donttest{stree=read.tree(text="((((a:10000,b:10000):10000,c:20000):10000,d:30000):10000,e:40000);")
#' pops=c(15000,25000,10000,1,1,1,1,1,12000)
#' gts=read.tree(file=system.file("extdata","genetreeSample",package="MSCsimtester"))
#' rootedTriple(stree,pops,gts,"a","b","c")
#' rootedTriple(stree,pops,gts,"a","b","c",1000) }
#' @export
#'
rootedTriple = function(stree,
                        popSizes,
                        gtSample,
                        taxon1,
                        taxon2,
                        taxon3,
                        subsampleSize = FALSE)
{
  # limited checking of arguments
  if (sum(duplicated(c(taxon1, taxon2, taxon3))) == 1) {
    warning('Taxa must be different.')
    return(invisible())
  }

  if (class(popSizes) != "numeric") {
    warning("Population sizes must be numeric. Exiting.")
    return(invisible())
  }
  else {
    invPopSizes = lapply(popSizes, invertPopSize)
  }

  if (subsampleSize == FALSE)
    subsampleSize = length(gtSample)

  droppedtips = stree$tip.label[!stree$tip.label %in% c(taxon1, taxon2, taxon3)]
  striplet = drop.tip(stree, droppedtips, root.edge = 0)

  numTests = floor(length(gtSample) / subsampleSize)
  pv = vector(mode = "numeric", length = numTests)
  intbrnch = vector(mode = "numeric", length = numTests)

  t1 = read.tree(text = paste0("((", taxon1, ",", taxon2, ")", ",", taxon3, ");"))
  t2 = read.tree(text = paste0("((", taxon1, ",", taxon3, ")", ",", taxon2, ");"))
  t3 = read.tree(text = paste0("((", taxon3, ",", taxon2, ")", ",", taxon1, ");"))
  ttops = c(
    paste0("((", taxon1, ",", taxon2, ")", ",", taxon3, ");"),
    paste0("((", taxon1, ",", taxon3, ")", ",", taxon2, ");"),
    paste0("((", taxon3, ",", taxon2, ")", ",", taxon1, ");")
  )
  if (all.equal(striplet, t1, use.edge.length = FALSE))
  {
    abo = edgesabove(stree, taxon1, taxon2)
    abo2 = edgesabove(stree, taxon3, taxon2)
    abo = abo[!abo %in% abo2]
    cu = 0

    for (i in 1:length(abo))
    {
      cu = stree$edge.length[abo[i]] / popSizes[abo[i]] + cu
    }
    trip = c(2, 3)
    tita = 1
  }
  if (all.equal(striplet, t2, use.edge.length = FALSE))
  {
    abo = edgesabove(stree, taxon1, taxon3)
    abo2 = edgesabove(stree, taxon3, taxon2)
    abo = abo[!abo %in% abo2]
    cu = 0

    for (i in 1:length(abo))
    {
      cu = stree$edge.length[abo[i]] / popSizes[abo[i]] + cu
    }
    trip = c(1, 3)
    tita = 2
  }
  if (all.equal(striplet, t3, use.edge.length = FALSE)) {
    abo = edgesabove(stree, taxon3, taxon2)
    abo2 = edgesabove(stree, taxon1, taxon2)
    abo = abo[!abo %in% abo2]
    cu = 0
    # Coalescent units
    for (i in 1:length(abo))
    {
      cu = stree$edge.length[abo[i]] / popSizes[abo[i]] + cu
    }
    trip = c(2, 1)
    tita = 3
  }
  ntop = c(1, 2, 3)[-tita]
  for (h in 1:numTests)
  {
    Q = matrix(0, subsampleSize, 4) # array for 3 counts
    for (z in 1:subsampleSize)
    {
      t <- gtSample[[(h - 1) * subsampleSize + z]]
      t = drop.tip(t, droppedtips, root.edge = 0)

      if (all.equal(t, t1, use.edge.length = FALSE))
      {
        Q[z, 4] = 1
      }

      if (all.equal(t, t2, use.edge.length = FALSE))
      {
        Q[z, 4] = 2
      }
      if (all.equal(t, t3, use.edge.length = FALSE))
      {
        Q[z, 4] = 3
      }
    }
    s = Q
    s4 = s[, 4]
    c = rep(0, 3)
    c[1] = sum(s4 == 1)
    c[2] = sum(s4 == 2)
    c[3] = sum(s4 == 3)
    cc = exp(-cu) / 3
    d = -log(1.5 * (c[trip[1]] + c[trip[2]]) / subsampleSize)
    intbrnch[h] = d
    theo = c(subsampleSize * (1 - 2 * cc),
             subsampleSize * cc,
             subsampleSize * cc)
    sim = c(c[tita], c[ntop[1]], c[ntop[2]])
    pv[h] = pchisq(sum((theo - sim) ^ 2 / theo), 2, lower.tail = F)

  }
  if (numTests == 1)
  {
    obj = list()
    smmry <-
      matrix(
        c(
          subsampleSize * (1 - 2 * cc),
          subsampleSize * cc,
          subsampleSize * cc,
          c[tita],
          c[ntop[1]],
          c[ntop[2]]
        ),
        ncol = 3,
        byrow = TRUE
      )
    colnames(smmry) <- c(ttops[tita], ttops[ntop[1]], ttops[ntop[2]])
    rownames(smmry) <- c('Expected counts:', 'Sample counts:')
    smmry <- as.table(smmry)
    obj$TripletCounts = smmry
    obj$pv = pv[1]
    smmry2 <- matrix(c(cu, d), ncol = 1, byrow = TRUE)
    rownames(smmry2) <-
      c('Species tree internal edge length (cu):',
        'Estimated internal edge length (cu):')
    colnames(smmry2) <- c(' ')
    smmry2 <- as.table(smmry2)
    obj$InternalEdge = smmry2
    class(obj) = "rootedTripleOutput"

    return(obj)
  } else {
    str = sprintf(
      "Chi^2 test results on samples of size %d,\n taxa: %s, %s, %s",
      subsampleSize,
      taxon1,
      taxon2,
      taxon3
    )
    hist(
      pv,
      breaks = numTests,
      xlab = "p-values for rooted triple tests",
      xlim = c(0, 1),
      freq = TRUE,
      main = str
    )
    str = sprintf("Internal branch length estimates,\n taxa: %s, %s, %s",
                  taxon1,
                  taxon2,
                  taxon3)
    hh = hist(
      intbrnch,
      breaks = numTests / 2,
      main = str,
      xlab = "Coalescent units",
      freq = TRUE
    )
    xV = rep(cu, 2)
    yV = c(-.2, .05 * max(hh$density))
    lines(xV, yV, lwd = 4, col = "blue")
  }
}





# This function determines the edges above two taxa.
# It returns list of edge numbers, from the MRCA to the root
#
edgesabove = function(stree, tx1, tx2)
{
  MRCA = mrca(stree, full = FALSE) #compute full array so also works for tx1=tx2 (getMRCA requires distinct taxa)
  MRCA = MRCA[tx1, tx2]
  ancestry = MRCA
  flag = length(which(stree$edge[, 2] == MRCA))
  p = 1
  while (flag == 1)
  {
    ancestry[p + 1] = stree$edge[which(stree$edge[, 2] == ancestry[p])]
    flag = length(stree$edge[which(stree$edge[, 2] == ancestry[p + 1])])
    p = p + 1
  }
  edge.ancestry = 1
  for (i in  1:length(ancestry))
  {
    if (length(which(stree$edge[, 2] == ancestry[i])) == 1)
      edge.ancestry[i] =  which(stree$edge[, 2] == ancestry[i])
    else {
      edge.ancestry[i] = nrow(stree$edge) + 1
    }
  }
  return(edge.ancestry)
}


# invert a numerical population size and make it a function
invertPopSize = function(pp)
{
  xx <- local({
    function(x) {
      1 / (pp + x * 0)
    }
  })
  return(xx)
}



#' Print function for objects of class \code{ADtestOutput}.

#' @param x An object of class \code{ADtestOutput}, as produced by \code{ADtest} function.
#' @param ... Further arguments to be passed to or from other methods.
#'
#' @return NONE

#' @export
print.ADtestOutput <- function(x, ...) {
  if (class(x$Adtest) == "kSamples") {
    print(x$Adtest)
  }
  else{
    cat(sprintf("  Anderson-Darling tests performed: %g \n", x$Adtest))
  }
}


#' Print function for objects of class \code{rootedTripleOutput}.
#'
#' @details Print function for objects of class \code{rootedTripleOutput}.
#'
#' @param x An object of class \code{rootedTripleOutput}, as produced by \code{rootedTriple} function.
#' @param ... Further arguments to be passed to or from other methods.
#'
#' @return NONE
#'
#' @export
print.rootedTripleOutput <- function(x, ...) {
  cat("  Rooted triple topology counts on gene trees \n \n")
  print(x$TripletCounts)
  cat(sprintf("  Chisq p-value: %g \n", x$pv))
  cat(" \n  Branch length on rooted triple in species tree")
  print(x$InternalEdge)
}

#' Simulated gene tree dataset.
#'
#' A dataset of 10,000 gene trees on 5 taxa simulated under the MSC on a species tree.
#'
#' @details This simulated dataset was produced by SimPhy \insertCite{simphy}{MSCsimtester},
#' using the species tree
#'
#' ((((a:10000,b:10000):10000,c:20000):10000,d:30000):10000,e:40000);
#'
#' with population sizes
#'
#' c(15000,25000,10000,1,1,1,1,1,12000)
#'
#' and edges ordered by the \code{ape} function \code{read.tree}.
#'
#' @docType data
#'
#' @name genetreeSample
#'
#' @details File is accessed as  \code{system.file("extdata","genetreeSample",package="MSCsimtester")}, for example
#' using the ape command:
#'
#' \code{gts=read.tree(file=system.file("extdata","genetreeSample",package="MSCsimtester") )}
#'
#' @format A text file with 10,000 metric Newick gene trees on the taxa a,b,c,d,e
#'
#' @references
#' \insertRef{simphy}{MSCsimtester}
#'
NULL

