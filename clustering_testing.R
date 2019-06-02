#'install.packages('cluster')
#'install.packages('fpc')
#'install.packages('NbClust')
#'install.packages('factoextra')

library(factoextra)
library(NbClust)
library(cluster)
library(fpc)

wdbc = read.table(file = "wdbc.data", header = FALSE, 
                  sep=',', row.names=1, na.strings="?")
colnames(wdbc) = c("diagnosis",
                   paste(c(rep("mean",10), rep("SE",10), rep("worst",10)),
                         rep(c("radius", "texture", "perimeter", "area",
                               "smoothness", "compactness", "concavity",
                               "concave_points", "symmetry", "fractal_dimension"),3),
                         sep="_"))


set.seed(123)

data = as.data.frame(scale(wdbc[-1],T,T))
cls.vec = as.numeric(wdbc$diagnosis)

#plot the dataset with reduced dimensions
X11()
fviz_pca_ind(prcomp(data), title = "PCA - WDBC data", ellipse.type = "convex",
             habillage = wdbc$diagnosis, palette = "Dark2",
             geom = "point", ggtheme = theme_minimal(),
             legend = "bottom"
             )

#' Checking cluster tendency
#' although the data is classified into two groups, it's still good to check it
#' hopkins stat is far below .5 -> data is sig. clusterable, the graph also shows a cluster structure of the data
cls.tend = get_clust_tendency(data, n = nrow(data)-1, graph = T, gradient = list(low = "darkgreen",  mid = "white", high = "darkred"))
cls.tend$hopkins_stat
cls.tend$plot+labs(title="WDBC data")

#' we could also generate random data from the WDBC dataset for comparison
random.data <- apply(data, 2, 
                     function(x){runif(length(x), min(x), (max(x)))})
random.data <- scale(as.data.frame(random.data))

cls.tend.random = get_clust_tendency(random.data, n = nrow(data)-1, graph = T, gradient = list(low = "darkgreen",  mid = "white", high = "darkred"))
cls.tend.random$hopkins_stat # H = ~0.5 -> uniform distribution
cls.tend.random$plot+labs(title="Random data") # no cluster structure


## Optimal nr of clusters - in accordance with our knowledge of the data structure
# for hierarchical clustering
nbcls = NbClust(data, distance = "euclidean",
                min.nc = 2, max.nc = 10, 
                method = "ward.D", index ="all") 
fviz_nbclust(nbcls, ggtheme = theme_minimal()) 

## & for e.g. k-means and clara with different methods - generally 2 wins
# k-means
# Elbow method -> suggests 4 
fviz_nbclust(data, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method - k-means")
# Silhouette method -> suggests 2
fviz_nbclust(data, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method - k-means")
# Gap statistic -> suggests 2
fviz_nbclust(data, kmeans,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "Gap statistic method - k-means")

# clara
# Elbow method -> suggests 4 
fviz_nbclust(data, clara, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method - clara")
# Silhouette method -> suggests 2
fviz_nbclust(data, clara, method = "silhouette")+
  labs(subtitle = "Silhouette method - clara")
# Gap statistic -> suggests 2
fviz_nbclust(data, clara,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "Gap statistic method - clara")



#### CLUSTERING (using eclust)
## Non-hierarchical
## 1. Kmeans
km = eclust(data, "kmeans", k = 2, iter.max = 50, nstart = 150, algorithm = "Lloyd", graph = F)
table(km$cluster, wdbc$diagnosis)

X11()
fviz_cluster(km, ellipse.type = "convex",
             palette = "Dark2", geom = "point",
             ggtheme = theme_minimal(), main = "K-means cluster plot")

X11()
# silhouette plot and average sil widths within groups and overall
fviz_silhouette(km, palette = "Dark2", subtitle = "K-means clustering")
km.sil = km$silinfo$avg.width

# mean squared distance from the centers
mean(km$withinss)


## 2. PAM
pam = eclust(data, "pam", k = 2, graph = F)
table(pam$clustering, wdbc$diagnosis)

X11()
fviz_cluster(pam, ellipse.type = "convex",
             palette = "Dark", geom = "point",
             ggtheme = theme_minimal(), main = "PAM cluster plot")

X11()
# silhouette plot and average sil widths within groups and overall
fviz_silhouette(pam, palette = "Dark", subtitle = "PAM clustering")
pam.sil = pam$silinfo$avg.width

## 3. Clara
clara = eclust(data, "clara", k = 2, graph = F)
table(clara$clustering, wdbc$diagnosis)

X11()
fviz_cluster(clara, ellipse.type = "convex",
             palette = "Set2", geom = "point",
             ggtheme = theme_minimal(), main = "Clara cluster plot")

X11()
# silhouette plot and average sil widths within groups and overall
fviz_silhouette(clara, palette = "Set2", subtitle = "Clara clustering")
clara.sil = clara$silinfo$avg.width

### Hierarchical clustering
## 4. Agglomerative - hclust
#' code to see which hclust gives best results using corrected Rand index (external validation)
#' as silhouette (internal validation) seems to not be reliable in some cases, i.e. the widths are high, 
#' but the clustering is poor, e.g. for method "complete" and metric "euclidean"
#' definig minimum sizes of clusters would be a solution(?) -> haven't found a way to do that for h.clustering

# example:
'
hclust.ex = eclust(data, "hclust", k=2, graph = F, hc_metric="euclidean", hc_method="complete")
table(wdbc$diagnosis, hclust.ex$cluster) # apart from 2 obs each is in the same (1st) cluster
fviz_cluster(hclust.ex, ellipse.type = "convex", palette = "npg", geom = "point", ggtheme = theme_minimal(), 
             main = "Hclust example cluster plot")
fviz_silhouette(hclust.ex, palette = "npg", subtitle = "Hclust example")
hclust.ex$silinfo$avg.width
'
#' this could be done also for agnes/diana, but generally I've just used 
#' the same metric and method for the rest of h.clust algorithms
methods <- c("single", "complete", "average", "ward.D", "ward.D2")
metrics <- c("euclidean", "manhattan")

best.hc = NULL
#rand = rep(0,length(methods)*length(metrics))
#k = 0 

for(mth in methods){
  for(mtr in metrics){
    #k = k + 1
    tmp.hc = eclust(data, "hclust", k=2, graph = FALSE, hc_metric=mtr, hc_method=mth)
    tmp.sts = cluster.stats(dist(data),  cls.vec, tmp.hc$cluster)
    #rand[k] = tmp.sts$corrected.rand
    if(is.null(best.hc)){
      best.hc = tmp.hc
    }
    else{
      best.stats = cluster.stats(dist(data), cls.vec, best.hc$cluster)
      if(best.stats$corrected.rand < tmp.sts$corrected.rand){
        best.hc = tmp.hc
      }
    }
  }
}
rm(tmp.hc)
rm(tmp.sts)

best.hc #euclidean, ward.D
table(best.hc$cluster, wdbc$diagnosis)

X11()
fviz_dend(best.hc, k = 2, show_labels = FALSE, palette = "uchicago", rect = TRUE, main = "Hclust Cluster Dendrogram")

X11()
fviz_cluster(best.hc, ellipse.type = "convex",
             palette = "uchicago", geom = "point",
             ggtheme = theme_minimal(), main = "Hclust cluster plot")

X11()
# silhouette plot and average sil widths within groups and overall
fviz_silhouette(best.hc, palette = "uchicago", subtitle = "Hclust hierarchical clustering")
hc.sil = best.hc$silinfo$avg.width


## 5. Agglomerative - agnes
agnes = eclust(data, "agnes", k = 2, graph = F, hc_metric = "euclidean", hc_method = "ward.D")
table(agnes$cluster, wdbc$diagnosis)

X11()
fviz_dend(agnes, k = 2, show_labels = FALSE, palette = "simpsons", rect = TRUE, main = "Agnes Cluster Dendrogram")

X11()
fviz_cluster(agnes, ellipse.type = "convex",
             palette = "simpsons", geom = "point",
             ggtheme = theme_minimal(), main = "Agnes cluster plot")

X11()
# silhouette plot and average sil widths within groups and overall
fviz_silhouette(agnes, palette = "simpsons", subtitle = "Agnes hierarchical clustering")
agnes.sil = agnes$silinfo$avg.width


## 6. Divisive hierarchical clustering - diana
diana = eclust(data, "diana", k=2, graph = FALSE, hc_metric="euc", hc_method="ward.D")
table(diana$cluster, wdbc$diagnosis)

X11()
fviz_dend(diana, k = 2, show_labels = FALSE, palette = "ucscgb", main = "Diana Cluster Dendrogram")

X11()
fviz_cluster(diana, ellipse.type = "convex",
             palette = "ucscgb", geom = "point",
             ggtheme = theme_minimal(), main = "Diana cluster plot")

X11()
# silhouette plot and average sil widths within groups and overall
fviz_silhouette(diana, palette = "ucscgb", subtitle = "Diana hierarchical clustering")
diana.sil = diana$silinfo$avg.width


## Quality measures
all.met = c('k-means', 'PAM', 'clara', 'agnes', 'hclust', 'diana')
# Internal validity - avg sil values for every method
avgSil = matrix(c(km.sil,pam.sil,clara.sil,agnes.sil,hc.sil,diana.sil))
colnames(avgSil) = c('Avg.sil.width')
rownames(avgSil) = all.met
avgSil = as.table(avgSil)
avgSil # the highest average silhouette width has diana

# all stats for each method
km.sts = cluster.stats(dist(data),  cls.vec, km$cluster)
pam.sts = cluster.stats(dist(data),  cls.vec, pam$clustering)
clara.sts = cluster.stats(dist(data),  cls.vec, clara$clustering)
agnes.sts = cluster.stats(dist(data),  cls.vec, agnes$cluster)
hc.sts = cluster.stats(dist(data),  cls.vec, best.hc$cluster)
diana.sts = cluster.stats(dist(data),  cls.vec, diana$cluster)

# External validity - corrected rand and vi measures
rand.vi = matrix(c(km.sts$corrected.rand,pam.sts$corrected.rand,clara.sts$corrected.rand,
                   agnes.sts$corrected.rand,hc.sts$corrected.rand,diana.sts$corrected.rand,
                   km.sts$vi,pam.sts$vi,clara.sts$vi,agnes.sts$vi,hc.sts$vi,diana.sts$vi),ncol=2)
colnames(rand.vi) = c('Rand','vi')
rownames(rand.vi) = all.met
rand.vi = as.table(rand.vi)

quality.tab = cbind(avgSil,rand.vi)
quality.tab 
'> quality.tab
        Avg.sil.width      Rand        vi
k-means     0.3449740 0.6707206 0.5772294
PAM         0.3491337 0.6078971 0.6563673
clara       0.3485067 0.6587003 0.5784200
agnes       0.3393848 0.5750409 0.7004186 
hclust      0.3036489 0.6373912 0.6435125
diana       0.3624915 0.4862299 0.6458926'

#' As we know the class of each observation, I'd first look at the external cluster validation measures;
#' the k-means with high Rand and small VI indexes seems to be the best choice for this particular type
#' of data (and there is no problem of choosing the appropriate cluster number), the k-means is followed by clara; 
#' both methods perform better than hierarchical clustering methods of which the hclust has the best results; 
#' visual re-comparison of the WDBC data and different clustering results plots supports the choice:
X11()
fviz_pca_ind(prcomp(data), title = "WDBC data", ellipse.type = "convex",
             habillage = wdbc$diagnosis, palette = "Dark2",
             geom = "point", ggtheme = theme_minimal())
X11()
fviz_cluster(km, ellipse.type = "convex",
             palette = "Dark2", geom = "point",
             ggtheme = theme_minimal(), main = "K-means cluster plot")
X11()
fviz_cluster(clara, ellipse.type = "convex",
             palette = "Set2", geom = "point",
             ggtheme = theme_minimal(), main = "Clara cluster plot")
X11()
fviz_cluster(best.hc, ellipse.type = "convex",
             palette = "lancet", geom = "point",
             ggtheme = theme_minimal(), main = "Hclust cluster plot")


# we can also compute the purity measure on the clustering results to further test them
ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

ClusterPurity(km$cluster, cls.vec) # the highest purity
ClusterPurity(clara$clustering, cls.vec)
ClusterPurity(best.hc$cluster, cls.vec)
# each of the rest has worse and worse scores
ClusterPurity(pam$clustering, cls.vec)
ClusterPurity(agnes$cluster, cls.vec)
ClusterPurity(diana$cluster, cls.vec)


### Additional (internal) comparison of the clustering methods
library(clValid)

# best three chosen on the basis of external validation measures
clmethods <- c("hierarchical","kmeans","clara")
valid <- c("internal", "stability")

intern <- clValid(data, nClust = 2, 
                            clMethods = clmethods, validation = "internal", metric = "euclidean", method = "ward")
stab <- clValid(data, nClust = 2, 
                  clMethods = clmethods, validation = "stability", metric = "euclidean", method = "ward")

# k-means wins only in terms of connectivity
summary(intern)
'Optimal Scores:

             Score   Method       Clusters
Connectivity 66.2083 kmeans       2       
Dunn          0.0660 hierarchical 2       
Silhouette    0.3485 clara        2'

#' however it outperforms clara and hclust in three stability (consistency of clustering) measures:
#' average proportion of non-overlap (APN);
#' average distance (AD);
#' average distance between means (ADM)
summary(stab)
'Optimal Scores:

Score  Method Clusters
APN 0.0206 kmeans 2       
AD  5.7918 kmeans 2       
ADM 0.1618 kmeans 2       
FOM 0.8259 clara  2 '


# and comparison of them all
clmethods <- c("hierarchical","kmeans","pam","diana","clara","agnes")
intern.all <- clValid(data, nClust = 2, 
                  clMethods = clmethods, validation = "internal", metric = "euclidean", method = "ward")
stab.all <- clValid(data, nClust = 2, 
                clMethods = clmethods, validation = "stability", metric = "euclidean", method = "ward")

# k-means has the best connectivity, agnes Dunn index and diana silhouette 
summary(intern.all)
'Optimal Scores:

Score   Method Clusters
Connectivity 66.2083 kmeans 2       
Dunn          0.0724 agnes  2       
Silhouette    0.3625 diana  2 '

# for stability measures, k-means still has the best result in APN and AD
summary(stab.all)
'Optimal Scores:

    Score  Method Clusters
APN 0.0206 kmeans 2       
AD  5.7918 kmeans 2       
ADM 0.1471 pam    2       
FOM 0.8259 clara  2 '
