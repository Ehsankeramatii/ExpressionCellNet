
############################################################################
############################################################################
############################################################################


# Start

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    crayon::green("ExpressionCellNet v1.0 loaded, written by Ehsan Keramati")
  )
}

.onLoad <- function(libname, pkgname) {
  needed <- c("igraph", "ggplot2", "crayon", "Metrics", "uwot")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "),
      "\nInstall them with:\n  install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "),
      "))",
      call. = FALSE
    )
  }
}

# Setup Functions :
############################################################################

########## Make ExpressionCellNEy Object ##########
#' Create an ExpressionCellNet object
#'
#' Wraps the core inputs (count matrix, annotation table, and interaction database)
#' into a single ExpressionCellNet object used by downstream functions.
#'
#' @param count_matrix A gene-by-sample expression count matrix (rows: genes, cols: samples).
#' @param annotations A data.frame containing gene annotations/mapping information.
#' @param intraction_db A data.frame of gene-gene interactions (two columns).
#'
#' @return An ExpressionCellNet object (list) containing core data components.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- createExpCellNetObj(
#'   count_matrix  = counts,
#'   annotations   = annot_df,
#'   intraction_db = interactions
#' )
#' }
createExpCellNetObj <- function(count_matrix , annotations , intraction_db) {

  # Train-Test Data Split :

  train_number <- round(ncol(count_matrix) * 0.8 , 0)
  test_number <- ncol(count_matrix) - train_number

  train_samples <- sample(1:ncol(count_matrix) , train_number , replace = F)
  test_samples <- setdiff(1:ncol(count_matrix) , train_samples )

  train_df <- data.frame(col1 = colnames(count_matrix)[train_samples] , col2 = train_samples , col3 = rep("Train" , train_number))
  test_df <- data.frame(col1 = colnames(count_matrix)[test_samples] , col2 = test_samples , col3 = rep("Test" , test_number))

  Train_test_split <- rbind(train_df , test_df)

  colnames(Train_test_split) <- c("Sample_Name" , "Sample_Number" , "Group")

  intcore_obj <- list(count_matrix ,annotations , intraction_db , Train_test_split)

  names(intcore_obj)[1] <- "Count Matrix"
  names(intcore_obj)[2] <- "Annotations"
  names(intcore_obj)[3] <- "Interaction DataBase"
  names(intcore_obj)[4] <- "Train-Test Samples"

  return(intcore_obj)

}

########## Count Matrix Normalization ##########
#' Normalize expression count matrix
#'
#' Filters low-expression genes (optional), applies CPM normalization (optional),
#' and log2 transformation (optional). Stores results inside the ExpressionCellNet object.
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#' @param CheckExpression Logical; if TRUE, low-expression genes are filtered.
#' @param CPM Logical; if TRUE, counts are converted to counts-per-million.
#' @param LogTransform Logical; if TRUE, log2 transform is applied (typically after CPM).
#'
#' @return Updated ExpressionCellNet object containing normalized matrix.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- CountMatrixNormalization(obj, CheckExpression = TRUE, CPM = TRUE, LogTransform = TRUE)
#' }
#'
CountMatrixNormalization <- function(ExpCellNetObj , CheckExpression = TRUE , CPM = TRUE , LogTransform = FALSE) {

  # Check Expression :

  df <- ExpCellNetObj[["Count Matrix"]]
  gene_number_before <- nrow(df)


  if (CheckExpression == TRUE) {

    cat(crayon::red("Removing genes that are low or not expressed...\n"))

    #n <- ncol(df) / 2
    #df <- df[which(rowSums(df) > n) ,]


    min_samples <- ncol(df) / 2
    df <- df[rowSums(df > 0) >= min_samples, ]



    cat(crayon::green("Done!\n"))

  }

  gene_number_after <- nrow(df)
  number_of_gene_deleted <- gene_number_before - gene_number_after
  cat(crayon::blue(paste0(number_of_gene_deleted , " genes with no detectable expression were excluded.\n")))

  # CPM :

  if (CPM == TRUE) {


    cat(crayon::red("Normalizing count matrix by CPM method...\n"))
    df <- apply(df, 2, function(x) (x / sum(x)) * 1e6)
    cat(crayon::green("Done!\n"))

  }

  # log2 Transform :

  if (LogTransform == TRUE) {


    cat(crayon::red("Transforming Data to Log2...\n"))
    df <- log2(df + 1)
    cat(crayon::green("Done!\n"))


  }

  df <- as.data.frame(df)

  if (length(ExpCellNetObj[["Normalized Count Matrix"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(df))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Normalized Count Matrix"
    return(ExpCellNetObj)
  }

  if (length(ExpCellNetObj[["Normalized Count Matrix"]]) != 0){

    ExpCellNetObj[["Normalized Count Matrix"]] <- df

    return(ExpCellNetObj)
  }



}

########## Bar plot for Count Matrix Normalization ##########

#' Plot expression matrix distribution
#'
#' Visualizes the distribution of raw or normalized expression values across samples.
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#' @param Matrix Character; which matrix to plot (e.g., "RAW" or "Normalized").
#'
#' @return A plot (base R) is produced as a side effect.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' PlotCountMatrix(obj, Matrix = "RAW")
#' PlotCountMatrix(obj, Matrix = "Normalized")
#' }
#'
PlotCountMatrix <- function(ExpCellNetObj , Matrix = "Normalized") {

  # Default Values :
  if (Matrix == "") {Matrix = "Normalized" }


  if (Matrix == "RAW") {

    df <- ExpCellNetObj[["Count Matrix"]]
    plot1 <- graphics::barplot(colSums(df) , main = "Total counts for each sample in rAW matrix",
                     col = "darkblue" , xlab = "Samples" , ylab = "Count",space = 1,border = F,names.arg = F)


  }

  if (Matrix == "Normalized") {

    df <- ExpCellNetObj[["Normalized Count Matrix"]]
    plot1 <- graphics::barplot(colSums(df) , main = "Total counts for each sample in normalized matrix" ,
                     col = "maroon" , xlab = "Samples" , ylab = "Count",space = 1,border = F,names.arg = F)

  }




}

# Build Network Functions :
############################################################################
#' Construct a gene interaction network centered on a seed gene
#'
#' Builds a multi-generation network starting from a seed gene using curated interactions
#' and correlation filtering across generations.
#'
#' @param ExpCellNetObj An ExpressionCellNet object containing expression and interaction data.
#' @param Gene Character; seed gene name (generation 0).
#' @param Generation Integer; number of generations to expand from the seed gene.
#' @param CorrelationThreshold Numeric vector; per-generation correlation thresholds.
#'
#' @return Updated ExpressionCellNet object with added components such as "Network" and "Generations".
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- BuildNetwork(obj, Gene = "PYGO2", Generation = 3, CorrelationThreshold = c(0.7,0.7,0.7))
#' }
#'
########## Build Network ##########

BuildNetwork <- function(ExpCellNetObj, Gene = "PYGO2" , Generation = 3 , CorrelationThreshold = c(0.7 , 0.7 , 0.7)) {


  count_df <- ExpCellNetObj[["Normalized Count Matrix"]]
  annot_df <- ExpCellNetObj[["Annotations"]]
  gene_symbol <- Gene
  interaction_DB <- ExpCellNetObj[["Interaction DataBase"]]
  gen <- Generation
  cor_threshold <- CorrelationThreshold
  meta_df <- c()
  meta_df <- rbind(meta_df , c(Gene , "Generation : 0"))



  fx <- function(count_df , annot_df , gene_symbol, interaction_DB ,cor_threshold) {


    string_list <- interaction_DB
    gene_ENS <- c()
    my_gene_network <- c()
    my_gene_network <-c()
    vec1 <- c()
    vec2 <- c()
    GOI_ID <- c()
    GOI_row_num <- c()


    gene_ENS <- annot_df[annot_df$SYMBOL == gene_symbol,3]
    my_gene_network <- string_list[which(string_list$V1 == gene_ENS | string_list$V2 == gene_ENS),]
    my_gene_network <- unique(my_gene_network)

    vec1 <- append(my_gene_network$V1 , my_gene_network$V2)
    vec1 <- unique(vec1)

    vec2 <- annot_df[annot_df$ENSEMBL %in% vec1 ,1]


    count_df <- count_df[rownames(count_df) %in% vec2 ,]


    GOI_ID <- annot_df[annot_df$SYMBOL ==gene_symbol , 1]
    GOI_row_num <- which(rownames(count_df) == GOI_ID)

    vec3 <- c()
    for( i in 1:nrow(count_df)) {

      new_cor <- stats::cor.test( as.numeric(count_df[i,]) , as.numeric(count_df[GOI_row_num,]) , method = "spearman")

      df <- data.frame(col1 = as.numeric(count_df[i,]) , col2 = as.numeric(count_df[GOI_row_num,]))

      # ---> Edited : Linear Model Was Removed

      #model <- lm(col1 ~ col2, data=df)
      #intercept <- coef(model)[1]
      #slope <- coef(model)[2]

      # ---> Edited : intercept and slope was removed
      vec3 <- rbind(vec3 , c(rownames(count_df)[i]  , rownames(count_df)[GOI_row_num] ,
                             new_cor$statistic , new_cor$estimate , new_cor$p.value))




    }

    # ---> Edited : intercept and slope was removed
    colnames(vec3) <- c("gene1" , "gene2" , "statistic" , "rho" , "pvalue")
    vec3 <- as.data.frame(vec3)
    #View(vec3)
    #hist(vec3$rho)

    # ---> Edited : intercept and slope was removed
    vec3$statistic <- as.numeric(vec3$statistic)
    vec3$rho <- as.numeric(vec3$rho)
    vec3$pvalue <- as.numeric(vec3$pvalue)

    vec4 <- subset(vec3 , vec3$pvalue < 0.01 & abs(vec3$rho) > cor_threshold)


    vec5 <-c()
    for( i in 1:nrow(vec4)) {

      a <- as.character(vec4[i,1])
      b <- as.character(vec4[i,2])
      vec5 <- rbind(vec5 , c(annot_df[annot_df$ID == a , 2] , annot_df[annot_df$ID == b , 2]  ))


    }


    # ---> Edited : before : 3:7 , Now : 3:5
    vec5 <- as.data.frame(vec5)
    vec5 <- cbind(vec5 , vec4[,c(3:5)])
    colnames(vec5)[1] <- "gene1"
    colnames(vec5)[2] <- "gene2"

    vec6 <-c()
    for(i in 1:nrow(vec5)) {

      i = 1
      a <- vec5[i,4]

      if(a > 0){

        vec6 <- append(vec6 , "+")
      }

      if(a < 0){

        vec6 <- append(vec6 , "-")
      }

    }

    # ---> Edited : before : 8 , now : 6
    vec5 <- cbind(vec5 , vec6)
    colnames(vec5)[6] <- "cor"

    #View(vec5)
    return(vec5)

  }


  F_0 <-c()
  genes_list <- c()
  genes_list <- c()

  print(paste0("Generaion : 1" , "/" , gen))
  F_0 <- fx(count_df , annot_df , gene_symbol,interaction_DB, cor_threshold[1])

  if ( nrow(F_0) < 2) {

    print("Cant detect any significant correlation in first generation")

  }


  if ( nrow(F_0) > 1) {


    genes_list <- F_0$gene1
    genes_list <- genes_list[genes_list != gene_symbol]

    print(paste0("Gene : " , gene_symbol ))

    F_j <- c()
    for(j in 2:gen) {


      print(paste0("Generaion : ", j , "/" , gen))



      data <- count_df  # note
      for(i in 1:length(genes_list)) {


        F_n <- fx(data , annot_df , genes_list[i],interaction_DB , cor_threshold[j])

        if (nrow(F_n) > 1 ) {

          F_j <- rbind(F_j , F_n)
          print(paste0("Gene : " , i , "/" , length(genes_list) , " - " , genes_list[i] ))
          meta_df <- rbind(meta_df ,c(genes_list[i] , paste0("Generation : ", j -1)) )
        }

        if (nrow(F_n) < 2 ) {


          print(paste0( "Gene : " , i , "/" , length(genes_list) , " - " , genes_list[i] , " in generation " , j , " does not have significant correlation with its neighbors"))
          meta_df <- rbind(meta_df ,c(genes_list[i],  paste0("Generation : ", j -1) ) )
        }

      }

      genes_list_2 <- F_j$gene1
      genes_list_2 <- unique(genes_list_2)
      genes_list_2 <- genes_list_2[genes_list_2 != gene_symbol]
      genes_list_2 <- genes_list_2[genes_list_2 %in% setdiff(genes_list_2 , genes_list)]
      genes_list <- genes_list_2


    }
    F_j <- rbind(F_j  , F_0)
    F_j <- subset(F_j , F_j$rho != 1)

    # remove duplicate edges

    rownames(F_j) <- NULL

    vec_x <- c()
    for(i in 1:nrow(F_j)) {

      a <- F_j[i,1]
      b <- F_j[i,2]

      c <- rownames(F_j[F_j$gene1 == b & F_j$gene2 == a,])

      if (length(c) == 1){

        c <- as.numeric(c)

        if ( c > i ) {
          vec_x <- append(vec_x , c)

        }

      }

      if (length(c) > 1){

        c <- c[1]
        c <- as.numeric(c)

        if ( c > i ) {
          vec_x <- append(vec_x , c)

        }

      }



    }

    vec_x <- as.numeric(vec_x)

    F_j <- F_j[-vec_x,]

    # unique :


    rownames(F_j) <- NULL

    vec_x <- c()

    for(i in 1:nrow(F_j)) {

      a <- F_j[i,1]
      b <- F_j[i,2]

      vec_x <- append(vec_x , rownames(F_j[F_j$gene1 == a & F_j$gene2 == b,])[1])

    }

    vec_x <- as.numeric(vec_x)
    F_j <- F_j[vec_x,]

    F_j <- unique(F_j)

    F_j <- as.data.frame(F_j)

    meta_df <- as.data.frame(meta_df)
    colnames(meta_df) <- c("Gene" , "Generation")

    all_gene <- c(F_j[,1] , F_j[,2])
    all_gene <- unique(all_gene)
    length(all_gene)

    genes_in_meta <- c(meta_df[,1])
    length(genes_in_meta)

    dif_genes <- setdiff(all_gene , genes_in_meta)

    new_df_2 <- data.frame(col1 = dif_genes , col2 = rep( paste0("Generation : ", Generation) , length(dif_genes)))
    colnames(new_df_2) <-  c("Gene" , "Generation")
    meta_df <- rbind(meta_df ,new_df_2)
    meta_df <- as.data.frame(meta_df)

    if(length(ExpCellNetObj[["Network"]]) == 0) {


      ExpCellNetObj <- append(ExpCellNetObj , list(F_j))
      names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Network"
      ExpCellNetObj <- append(ExpCellNetObj , list(meta_df))
      names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Generations"
      return(ExpCellNetObj)

    }

    if(length(ExpCellNetObj[["Network"]]) != 0) {


      ExpCellNetObj[["Network"]] <- F_j
      ExpCellNetObj[["Generations"]] <- meta_df

      return(ExpCellNetObj)

    }



  }

}

########## Show Network ##########
#' Visualize the constructed network
#'
#' Plots the current network stored in the ExpressionCellNet object.
#'
#' @param ExpCellNetObj An ExpressionCellNet object containing a "Network" edge table.
#'
#' @return A plot (igraph/base) is produced as a side effect.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- BuildNetwork(obj, Gene="PYGO2", Generation=3, CorrelationThreshold=c(0.7,0.7,0.7))
#' ShowNetwork(obj)
#' }
#'
ShowNetwork <- function(ExpCellNetObj){



  df <- ExpCellNetObj[["Network"]]
  new_graph <- igraph::graph_from_data_frame(df[,c(1,2)], directed = F )

  Generations <- ExpCellNetObj[["Generations"]]
  Main_node <- Generations[1,1]


  igraph::V(new_graph)$color <- "aquamarine2"
  igraph::V(new_graph)[igraph::V(new_graph)$name %in% Main_node]$color <- "hotpink1"
  plot(new_graph ,  vertex.label.cex = 0.6)

}

########## Find Pathway ##########

#' Find the shortest path between two genes in the network
#'
#' Computes the shortest interaction path (as a sequence of gene names) between a
#' start gene and an end gene using the network stored in the ExpressionCellNet object.
#'
#' @param ExpCellNetObj An ExpressionCellNet object containing a constructed network.
#' @param Start_gene Character; source gene symbol (must exist as a node in the network).
#' @param End_gene Character; destination gene symbol (must exist as a node in the network).
#'
#' @return A character vector of gene symbols representing the shortest path.
#'
#' @export
#' @importFrom igraph graph_from_data_frame shortest_paths V
#'
#' @examples
#' \dontrun{
#' obj <- BuildNetwork(obj, Gene="PYGO2", Generation=3, CorrelationThreshold=c(0.7,0.7,0.7))
#' path <- FindPathway(obj, Start_gene="PYGO2", End_gene="ATAD5")
#' path
#' }


FindPathway <- function(ExpCellNetObj, Start_gene , End_gene){

  df <- ExpCellNetObj[["Network"]]

  new_graph <- igraph::graph_from_data_frame(
    df[, 1:2],
    directed = FALSE
  )

  path <- igraph::shortest_paths(
    new_graph,
    from = Start_gene,
    to = End_gene
  )$vpath[[1]]

  path_names <- igraph::V(new_graph)[path]$name

  return(path_names)

}

########## Show Pathway ##########

#' Visualize a pathway (highlighted subgraph)
#'
#' Highlights the nodes/edges of a provided pathway on the full network visualization.
#'
#' @param ExpCellNetObj An ExpressionCellNet object containing a constructed network.
#' @param Pathway A character vector of genes representing a path (e.g., output of FindPathway()).
#'
#' @return A plot (igraph/base) is produced as a side effect.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' p <- FindPathway(obj, Start_gene="PYGO2", End_gene="ATAD5")
#' ShowPathway(obj, Pathway = p)
#' }
#'
ShowPathway <- function(ExpCellNetObj , Pathway) {


  df <- ExpCellNetObj[["Network"]]

  new_graph <- igraph::graph_from_data_frame(
    df[, 1:2],
    directed = FALSE
  )

  path_names <- Pathway
  nodes_to_color <- path_names

  igraph::V(new_graph)$color <- "aquamarine2"
  igraph::V(new_graph)[igraph::V(new_graph)$name %in% nodes_to_color]$color <- "tomato"

  igraph::plot.igraph(
    new_graph,
    vertex.label.cex = 0.6
  )


}

# Regression Prediction Functions : :
############################################################################


########## Plot Pathway Prediction ##########

#' Plot pathway-restricted prediction to a target gene
#'
#' Summarizes modeled expression changes along the seed-to-target pathway.
#'
#' @param ExpCellNetObj An ExpressionCellNet object after NetworkPrediction().
#' @param target_gene Character; destination gene to focus on.
#'
#' @return Invisibly returns NULL (plot is produced).
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar geom_text labs scale_fill_manual
#'
#' @examples
#' \dontrun{
#' obj <- NetworkPrediction(obj, Gene="PYGO2", Expression=0)
#' PlotPathwayPrediction(obj, "ATAD5")
#' }

PlotPathwayPrediction <- function(ExpCellNetObj , target_gene) {
  vec0 <- ExpCellNetObj[["Network Regression Analysis Results"]]

  First_gene <- vec0[1, 1]
  last_gene <- target_gene

  Pathway <- FindPathway(ExpCellNetObj, First_gene, last_gene)

  vec1 <- vec0[vec0[, 1] %in% Pathway, ]

  k <- c(
    "tomato", "maroon", "aquamarine4", "darkolivegreen", "purple1", "seagreen2",
    "khaki", "cyan4", "darkblue", "royalblue", "hotpink2", "purple4", "honeydew3",
    "steelblue1", "slateblue3", "lightsalmon"
  )

  x1 <- c()
  if (base::mean(base::as.numeric(vec1[, 4])) > 0) {
    x1 <- -1
  }
  if (base::mean(base::as.numeric(vec1[, 4])) < 0) {
    x1 <- 1
  }

  p1 <- ggplot2::ggplot(
    vec1,
    ggplot2::aes(
      x = vec1[, 1],
      y = base::as.numeric(vec1[, 4]),
      fill = vec1[, 1]
    )
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(
      title = "Gene Expression Changes",
      x = "Genes",
      y = "Changes (%)"
    ) +
    ggplot2::scale_fill_manual(values = base::sample(k, base::nrow(vec1))) +
    ggplot2::labs(fill = "Genes") +
    ggplot2::geom_text(
      ggplot2::aes(label = base::round(base::as.numeric(vec1[, 4]), 2)),
      vjust = x1,
      hjust = 1.2
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = base::as.numeric(vec1[,8]),
        ymax = base::as.numeric(vec1[,9])
      ),
      width = 0.2
    )

  return(p1)
}

########## Network Linear Regression Prediction ##########
#' Predict network-wide expression changes after perturbing a seed gene
#'
#' Fits linear regression models for seed-to-network genes and predicts shifts under a user-defined
#' perturbation scenario. Stores multiple result components in the object.
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#' @param Gene Character; seed gene to perturb.
#' @param Expression Numeric; imposed expression value under perturbation.
#'
#' @return Updated ExpressionCellNet object with regression models and results.
#' @export
#'
#' @importFrom stats lm as.formula coef
#'
#' @examples
#' \dontrun{
#' obj <- NetworkPrediction(obj, Gene="PYGO2", Expression=0)
#' }
NetworkPrediction <- function(ExpCellNetObj , Gene , Expression) {

  # Input Data :

  #test input
  #ExpCellNetObj <- intcoreobj
  #Gene <- "PYGO2"
  #Expression <- 0
  #

  df <- ExpCellNetObj[["Network"]]
  annot_df <- ExpCellNetObj[["Annotations"]]
  Count_matrix <- ExpCellNetObj[["Normalized Count Matrix"]]

  feature_df <- ExtractFeature(ExpCellNetObj)

  # Gene made nazar bayad biad avale jadval

  feature_df_sorted <- c()
  feature_df_sorted <- cbind(feature_df_sorted ,feature_df[,colnames(feature_df) == Gene] )
  colnames(feature_df_sorted)[1] <- Gene

  dif_genes <- setdiff(colnames(feature_df) ,Gene )
  feature_df <- feature_df[,colnames(feature_df) %in% dif_genes ]

  feature_df_sorted <- cbind(feature_df_sorted , feature_df)
  feature_df_sorted <- as.data.frame(feature_df_sorted)

  # Train #Test

  number_of_train_data <-round(nrow(feature_df_sorted) * 0.8 , 0)
  number_of_test_data <- nrow(feature_df_sorted) - number_of_train_data


  tarin_test_split_df <- ExpCellNetObj[["Train-Test Samples"]]

  train_data_rows <- tarin_test_split_df[tarin_test_split_df$Group == "Train" , 2]
  test_data_row <- tarin_test_split_df[tarin_test_split_df$Group == "Test" , 2]

  train_data <- feature_df_sorted[train_data_rows,]
  test_data <- feature_df_sorted[test_data_row,]


  # baraye inke majboor nasham kole codaye payin ro taghir bedam
  feature_df_sorted <- train_data


  models <- list()
  for (i in 2:ncol(feature_df_sorted)) {
    formula_str <- paste0(colnames(feature_df_sorted)[i], " ~ " , colnames(feature_df_sorted)[1])
    models[[i - 1]] <- stats::lm(stats::as.formula(formula_str), data = feature_df_sorted)
  }


  # add models to main object

  if (length(ExpCellNetObj[["Network Regression Models"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(models))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Network Regression Models"

  }

  if (length(ExpCellNetObj[["Network Regression Models"]]) != 0){

    ExpCellNetObj[["Network Regression Models"]] <- models

  }


  result_df <-c()
  result_df <- rbind(result_df , c(Gene ,Expression ))
  for(i in 1:length(models)) {

    slope = stats::coef(models[[i]])[2]
    intercept <- stats::coef(models[[i]])[1]

    X = Expression
    Y <- X*slope + intercept

    new_result <- c(colnames(feature_df_sorted)[i + 1] , Y )
    result_df <- rbind(result_df ,new_result )

  }

  rownames(result_df) <- NULL
  colnames(result_df) <- c("Gene" , "New_Expression")



  mean_exp <- c()
  for(i in 1:nrow(result_df)) {

    a <- result_df[i,1]
    b <- MeanExpression(ExpCellNetObj , a)
    mean_exp <- append(mean_exp , b)
  }


  result_df <- cbind(result_df , mean_exp)
  colnames(result_df)[3] <- "Mean_Expression"


  vec <- c()
  for(i in 1:nrow(result_df)){

    a <- as.numeric(result_df[i,2])
    b <- as.numeric(result_df[i,3])

    c <- (a-b)/b * 100

    if (c == "Inf"){ c = 100}
    if (c == "-Inf"){ c = -100}

    c <- round(c , 1)

    vec <- append(vec ,c )

  }

  result_df <- cbind(result_df , vec)
  colnames(result_df)[4] <- "Changes (%)"

  result_df_sorted <- result_df[,c(1,3)]
  result_df_sorted <- cbind(result_df_sorted ,result_df[,c(2,4)])


  result_df_sorted <- as.data.frame(result_df_sorted)


  # Testing Model


  test_result <-c()
  for(i in 1:length(models)) {

    #print(paste0("i = " , i))

    target_gene <- models[[i]][["terms"]][[2]]

    slope = stats::coef(models[[i]])[2]
    intercept <- stats::coef(models[[i]])[1]

    vec <- c()
    for(j in 1:nrow(test_data)) {


      #print(paste0("j = " , j))
      X = test_data[j,1]
      Y <- X*slope + intercept
      Y <- as.numeric(Y)
      vec <- append(vec ,Y )

    }

    test_result <- cbind(test_result , vec)
    test_result <- as.data.frame(test_result)

  }


  colnames(test_result) <- colnames(test_data)[2:ncol(test_data)]
  rownames(test_result) <- rownames(test_data)

  # test_data = Actual
  # test_result = pred

  real_data <- test_data[,-1]


  #View(real_data)
  #View(test_result)

  # add Actual Values to main object

  real_data <- as.data.frame(real_data)

  if (length(ExpCellNetObj[["Network Actual Values"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(real_data))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Network Actual Values"

  }

  if (length(ExpCellNetObj[["Network Actual Values"]]) != 0){

    ExpCellNetObj[["Network Actual Values"]] <- real_data

  }


  # add Predicted Values to main object

  test_result <- as.data.frame(test_result)

  if (length(ExpCellNetObj[["Network Predicted Values"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(test_result))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Network Predicted Values"

  }

  if (length(ExpCellNetObj[["Network Predicted Values"]]) != 0){

    ExpCellNetObj[["Network Predicted Values"]] <- test_result

  }


  accuracy_df <- c()
  vec <- c()
  RSE <- c()
  RMSE <- c()
  MAE <- c()
  MAPE <- c()
  RAE <- c()

  for( i in 1:ncol(test_result)) {

    pred <- test_result[,i]
    actual <- real_data[,i]


    gene_ens <- annot_df[annot_df[,2] == colnames(test_result)[i],1]
    gene_range <- Count_matrix[rownames(Count_matrix) == gene_ens,]

    min_data <- min(gene_range)
    max_data <- max(gene_range)

    data_scale <- max_data - min_data

    acc <- round(mean(100 - abs(pred -actual)/data_scale * 100),2)




    new_rse <- Metrics::rse(actual , pred)
    RSE <- append(RSE , new_rse)

    new_rmse <- Metrics::rmse(actual , pred)
    RMSE <- append(RMSE  , new_rmse)

    new_mae <- Metrics::mae(actual , pred)
    MAE <- append(MAE , new_mae)

    new_mape <- Metrics::mape(actual , pred)
    MAPE <- append(MAPE , new_mape)

    new_rae <- Metrics::rae(actual , pred)
    RAE <- append(RAE , new_rae)



    vec <- rbind(vec , c(colnames(test_result)[i] , acc))

  }


  RSE <- c(0 , RSE)
  RMSE <- c(0 , RMSE)
  MAE <- c(0 , MAE)
  MAPE <- c(0 , MAPE)
  RAE <- c(0 , RAE)


  colnames(vec) <- c("Gene" , "% Accuracy")

  mean_acc <- mean(as.numeric(vec[,2]))
  vec <- rbind(vec , c("Mean" , round(mean_acc , 2)))

  accuracy_df <- vec

  accuracy_df <- as.data.frame(accuracy_df)



  accuracy_df_2 <- accuracy_df[-nrow(accuracy_df),]


  result_df_sorted_2 <- result_df_sorted[-1,]

  result_df_sorted <- cbind(result_df_sorted ,c(0 ,accuracy_df[-nrow(accuracy_df),2] ) )

  colnames(result_df_sorted)[5] <- "Model Accuracy %"


  vec <- c()
  vec <-rbind(vec , c(0 , 0))
  for(i in 1:nrow(result_df_sorted_2)) {

    a <- as.numeric(result_df_sorted_2[i,3])
    b <- as.numeric(accuracy_df_2[i,2])
    c <- (100 - b) / 100
    d <- a*c

    a_min <- a - d
    a_max <- a + d

    vec <- rbind(vec , c(a_min ,a_max ))
  }

  colnames(vec) <- c("Min_Expression" , "Max_Expression")

  result_df_sorted <- cbind(result_df_sorted ,vec )

  vec <- c()
  vec <-rbind(vec , c(0 , 0))
  for(i in 1:nrow(result_df_sorted_2)) {

    a <- as.numeric(result_df_sorted_2[i,4])
    b <- as.numeric(accuracy_df_2[i,2])
    c <- (100 - b) / 100
    d <- a*c
    d <- abs(d)

    a_min <- a - d
    a_max <- a + d

    vec <- rbind(vec , c(a_min ,a_max ))
  }


  colnames(vec) <- c("Min_Percent" , "Max_Percent")
  result_df_sorted <- cbind(result_df_sorted ,vec )



  for(i in 1:nrow(result_df_sorted)) {

    a <- as.numeric(result_df_sorted[i,4])

    if (a > 100){ a = 100}
    if (a < -100){ a = -100}

    result_df_sorted[i,4] <- a

    b <- as.numeric(result_df_sorted[i,8])
    if (b > 100){ b = 100}
    if (b < -100){ b = -100}

    result_df_sorted[i,8] <- b


    c <- as.numeric(result_df_sorted[i,9])
    if (c > 100){ c = 100}
    if (c < -100){ c = -100}

    result_df_sorted[i,9] <- c

  }

  result_df_sorted[1,8] <- result_df_sorted[1,4]
  result_df_sorted[1,9] <- result_df_sorted[1,4]

  result_df_sorted <- as.data.frame(result_df_sorted)
  result_df_sorted[1,5] <- "-"

  if (length(ExpCellNetObj[["Network Regression Analysis Results"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(result_df_sorted))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Network Regression Analysis Results"

  }

  if (length(ExpCellNetObj[["Network Regression Analysis Results"]]) != 0){

    ExpCellNetObj[["Network Regression Analysis Results"]] <- result_df_sorted

  }



  Accuracy_df <- c()
  Accuracy_df <- cbind(result_df_sorted[,1] , result_df_sorted[,5] , RSE , RMSE , MAE , MAPE , RAE)

  colnames(Accuracy_df)[1] <- "Gene"
  colnames(Accuracy_df)[2] <- "Accuracy %"
  Accuracy_df[1,2:ncol(Accuracy_df)] <- "-"


  Accuracy_df <- Accuracy_df[-1,]
  Accuracy_df <- as.data.frame(Accuracy_df)


  if (length(ExpCellNetObj[["Network Analysis Accuracy Results"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(Accuracy_df))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Network Analysis Accuracy Results"

  }

  if (length(ExpCellNetObj[["Network Analysis Accuracy Results"]]) != 0){

    ExpCellNetObj[["Network Analysis Accuracy Results"]] <- Accuracy_df

  }



  return(ExpCellNetObj)

}

########## Plot Network Prediction ##########

#' Plot network prediction summary
#'
#' Creates a summary plot of predicted percent changes across network genes (bar plot).
#'
#' @param ExpCellNetObj An ExpressionCellNet object after NetworkPrediction().
#'
#' @return Invisibly returns NULL (plot is produced).
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar labs theme element_text scale_fill_gradient
#'
#' @examples
#' \dontrun{
#' obj <- NetworkPrediction(obj, Gene="PYGO2", Expression=0)
#' PlotNetworkPrediction(obj)
#' }

PlotNetworkPrediction <-function(ExpCellNetObj) {

  NetworkExpression <- ExpCellNetObj[["Network Regression Analysis Results"]]
  p1 <- ggplot2::ggplot(
    NetworkExpression,
    ggplot2::aes(
      x    = NetworkExpression[, 1],
      y    = as.numeric(NetworkExpression[, 4]),
      fill = as.numeric(NetworkExpression[, 4])
    )
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
    ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(
      title = "Gene Expression Changes",
      x     = "Genes",
      y     = "Changes (%)"
    ) +
    ggplot2::scale_fill_gradient(
      low  = "maroon",
      high = "khaki"
    ) +
    ggplot2::labs(fill = "Genes") +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = as.numeric(NetworkExpression[, 8]),
        ymax = as.numeric(NetworkExpression[, 9])
      ),
      width = 0.2
    )

  return(p1)
}


# Network Analysis Functions :
############################################################################

########## Network Analysis ##########

#' Analyze network topology and identify hub genes
#'
#' Computes degree, closeness, betweenness, and eigenvector centrality and stores:
#' (i) Network Analysis table and (ii) Hub Genes table (top 10%).
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#'
#' @return Updated ExpressionCellNet object with "Network Analysis" and "Hub Genes".
#' @export
#'
#' @importFrom igraph graph_from_data_frame degree closeness betweenness eigen_centrality
#'
#' @examples
#' \dontrun{
#' obj <- AnalyzeNetwork(obj)
#' }

AnalyzeNetwork <- function(ExpCellNetObj) {


  net_df <- ExpCellNetObj[["Network"]]

  g <- igraph::graph_from_data_frame(net_df[, 1:2], directed = FALSE)

  centrality_df <- data.frame(
    degree      = igraph::degree(g),
    closeness   = igraph::closeness(g),
    betweenness = igraph::betweenness(g),
    eigenvector = igraph::eigen_centrality(g)$vector
  )

  centrality_df <- centrality_df[order(centrality_df$degree, decreasing = TRUE), ]

  ExpCellNetObj[["Network Analysis"]] <- centrality_df

  ## Hub genes (top 10%)
  top_n <- ceiling(nrow(centrality_df) * 0.1)

  hubs <- Reduce(
    intersect,
    list(
      rownames(centrality_df)[order(centrality_df$degree, decreasing = TRUE)][1:top_n],
      rownames(centrality_df)[order(centrality_df$closeness, decreasing = TRUE)][1:top_n],
      rownames(centrality_df)[order(centrality_df$betweenness, decreasing = TRUE)][1:top_n],
      rownames(centrality_df)[order(centrality_df$eigenvector, decreasing = TRUE)][1:top_n]
    )
  )

  ExpCellNetObj[["Hub Genes"]] <- centrality_df[hubs, , drop = FALSE]

  return(ExpCellNetObj)

}

########## Show Hub Genes ##########

#' Show hub genes on the network
#'
#' Highlights hub genes on the igraph network visualization.
#'
#' @param ExpCellNetObj An ExpressionCellNet object after AnalyzeNetwork().
#'
#' @return Invisibly returns NULL (plot is produced).
#' @export
#'
#' @importFrom igraph graph_from_data_frame V
#'
#' @examples
#' \dontrun{
#' obj <- AnalyzeNetwork(obj)
#' ShowHubGenes(obj)
#' }

ShowHubGenes <- function(ExpCellNetObj) {

  df <- ExpCellNetObj[["Network"]]
  new_graph <- igraph::graph_from_data_frame(df[,c(1:2)], directed = F)
  hub_genes <- ExpCellNetObj[["Hub Genes"]]
  nodes_to_color <- rownames(hub_genes)
  igraph::V(new_graph)$color <- "aquamarine2"
  igraph::V(new_graph)[igraph::V(new_graph)$name %in% nodes_to_color]$color <- "khaki3"
  plot(new_graph ,  vertex.label.cex = 0.6)


}

# Enrichment Functions :
############################################################################

########## Gene Enrichment ##########

#' Functional enrichment analysis of network genes
#'
#' Runs enrichment analysis for genes in the constructed network using enrichR libraries.
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#'
#' @return A named list of enrichment result data.frames (per database).
#' @export
#'
#' @importFrom enrichR enrichr
#'
#' @examples
#' \dontrun{
#' enr <- GeneEnrichment(obj)
#' }

GeneEnrichment <- function(ExpCellNetObj){

  db <- c(
    "KEGG_2013", "KEGG_2015", "KEGG_2016", "KEGG_2019_Human", "KEGG_2021_Human",
    "Reactome_2013", "Reactome_2015", "Reactome_2016", "Reactome_2022", "Reactome_Pathways_2024",
    "WikiPathways_2013", "WikiPathways_2019_Human", "WikiPathways_2024_Human",
    "WikiPathways_2016", "Human_Gene_Atlas",
    "GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023",
    "GO_Biological_Process_2025", "GO_Cellular_Component_2025", "GO_Molecular_Function_2025"
  )

  df <- ExpCellNetObj[["Network"]]

  genes <- unique(c(df$gene1, df$gene2))

  enrichment <- enrichR::enrichr(genes, db)

  return(enrichment)

}

########## Bar plot for Enrichment ##########
#' Plot enrichment results
#'
#' Visualizes enrichment results for a selected database table.
#'
#' @param GeneEnrichment_dataframe A single enrichment result data.frame.
#' @param title Character; plot title.
#' @param order Character; ordering criterion (e.g., "pvalue").
#'
#' @return Invisibly returns NULL (plot is produced).
#' @export
#'
#' @importFrom enrichR plotEnrich
#'
#' @examples
#' \dontrun{
#' enr <- GeneEnrichment(obj)
#' PlotEnrichment(enr$KEGG_2021_Human, title="KEGG_2021_Human", order="pvalue")
#' }
PlotEnrichment <- function(GeneEnrichment_dataframe , title , order) {

  enrichR::plotEnrich(GeneEnrichment_dataframe, title = title , orderBy = order)

}

# Dimension Reduction Functions :
############################################################################

########## PCA for Network ##########

#' PCA on network-derived features
#'
#' Performs principal component analysis (PCA) on the network-restricted feature matrix
#' returned by ExtractFeature(), after min-max scaling, and returns a ggplot of PC1 vs PC2.
#'
#' @param ExpCellNetObj An ExpressionCellNet object containing "Network", "Generations",
#' and a normalized count matrix.
#'
#' @return A ggplot object showing PCA embedding of genes (nodes).
#'
#' @export
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_point geom_text labs theme
#'
#' @examples
#' \dontrun{
#' obj <- BuildNetwork(obj, Gene="PYGO2", Generation=3, CorrelationThreshold=c(0.7,0.7,0.7))
#' p <- PCAforGenes(obj)
#' p
#' }
#'
PCAforGenes <- function(ExpCellNetObj) {

  df <- ExpCellNetObj[["Network"]]
  annot_df <- ExpCellNetObj[["Annotations"]]
  Count_matrix <- ExpCellNetObj[["Normalized Count Matrix"]]

  all_features <- ExtractFeature(ExpCellNetObj)

  maxs <- base::apply(all_features, 2, base::max, na.rm = TRUE)
  mins <- base::apply(all_features, 2, base::min, na.rm = TRUE)

  scale_data <- function(x) {
    base::as.data.frame(
      base::scale(x, center = mins, scale = (maxs - mins))
    )
  }

  all_features <- scale_data(all_features)

  meta_pca <- ExpCellNetObj[["Generations"]]

  vec <- base::character(0)
  for (i in base::seq_len(base::ncol(all_features))) {
    a <- base::colnames(all_features)[i]
    b <- meta_pca[meta_pca[, 1] == a, 2][1]
    vec <- base::c(vec, b)
  }

  pc <- stats::prcomp(all_features)
  pcr <- base::data.frame(pc$rotation[, 1:3, drop = FALSE], vec, stringsAsFactors = FALSE)

  Groupe <- vec

  pca_plot <- ggplot2::ggplot(
    pcr,
    ggplot2::aes(
      x = pcr$PC1,
      y = pcr$PC2,
      color = Groupe,
      label = base::rownames(pcr)
    )
  ) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_text(size = 2, nudge_y = -0.01)

  pca_plot

}

########## UMAP for Network ##########
#' UMAP on network-derived features
#'
#' Computes a 2D UMAP embedding of the network-restricted feature matrix returned by
#' ExtractFeature() (after min-max scaling) and returns a ggplot visualization.
#'
#' @param ExpCellNetObj An ExpressionCellNet object containing "Network", "Generations",
#' and a normalized count matrix.
#'
#' @return A ggplot object showing UMAP embedding of genes (nodes).
#'
#' @export
#' @importFrom uwot umap
#' @importFrom ggplot2 ggplot aes geom_point geom_text labs theme element_blank
#'
#' @examples
#' \dontrun{
#' obj <- BuildNetwork(obj, Gene="PYGO2", Generation=3, CorrelationThreshold=c(0.7,0.7,0.7))
#' u <- UMAPforGenes(obj)
#' u
#' }
UMAPforGenes <- function(ExpCellNetObj) {

  df           <- ExpCellNetObj[["Network"]]
  annot_df     <- ExpCellNetObj[["Annotations"]]
  Count_matrix <- ExpCellNetObj[["Normalized Count Matrix"]]
  Generations  <- ExpCellNetObj[["Generations"]]

  all_features <- ExtractFeature(ExpCellNetObj)

  maxs <- base::apply(all_features, 2, base::max)
  mins <- base::apply(all_features, 2, base::min)

  scale_data <- function(df) {
    base::as.data.frame(
      base::scale(df, center = mins, scale = maxs - mins)
    )
  }

  all_features <- scale_data(all_features)
  all_features <- base::t(all_features)

  umap_result <- uwot::umap(all_features)

  vec <- c()
  for (i in seq_len(nrow(umap_result))) {
    a <- base::rownames(umap_result)[i]
    b <- Generations[Generations$Gene == a, 2]
    vec <- base::append(vec, b)
  }

  umap_result <- base::cbind(umap_result, vec)
  umap_result <- base::as.data.frame(umap_result)

  ggplot2::ggplot(
    umap_result,
    ggplot2::aes(
      x = umap_result[, 1],
      y = umap_result[, 2],
      color = umap_result[, 3],
      label = rownames(umap_result)
    )
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      title = "UMAP for Nodes",
      x = "UMAP2",
      y = "UMAP1",
      color = "Legendary"
    ) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(size = 2, nudge_y = 0.4)

}

# Add Data Functions :
############################################################################

########## Add Clinical Data ##########

#' Add clinical data
#'
#' Optionally attaches clinical variables to the ExpressionCellNet object.
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#' @param ClinicalData A data.frame of clinical variables (rows=samples).
#'
#' @return Updated ExpressionCellNet object.
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- AddClinicalData(obj, clinical_df)
#' }

AddClinicalData <- function(ExpCellNetObj , ClinicalData) {

  if (length(ExpCellNetObj[["Clinical Data"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(ClinicalData))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Clinical Data"
    return(ExpCellNetObj)
  }

  if (length(ExpCellNetObj[["Clinical Data"]]) != 0){

    ExpCellNetObj[["Clinical Data"]] <- list(ClinicalData)

    return(ExpCellNetObj)
  }

}

########## Add MetaData ##########

#' Add sample-level metadata
#'
#' Optionally attaches sample-level metadata to the ExpressionCellNet object.
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#' @param MetaData A data.frame of metadata (rows=samples).
#'
#' @return Updated ExpressionCellNet object.
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- AddMetaData(obj, meta_df)
#' }

AddMetaData <- function(ExpCellNetObj , MetaData) {

  if (length(ExpCellNetObj[["MetaData"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(MetaData))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "MetaData"
    return(ExpCellNetObj)
  }

  if (length(ExpCellNetObj[["MetaData"]]) != 0){

    ExpCellNetObj[["MetaData"]] <- list(MetaData)

    return(ExpCellNetObj)
  }

}

# Tools :
############################################################################

########## Mean Expression ##########

#' Compute mean expression of a gene
#'
#' Computes mean expression from raw or normalized matrix (default: normalized).
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#' @param Gene Character; target gene.
#' @param Normalized Logical; use normalized matrix (TRUE) or raw (FALSE).
#'
#' @return Numeric mean expression value.
#' @export

MeanExpression <- function(ExpCellNetObj , Gene , Normalized = TRUE ) {

  annot_df <- ExpCellNetObj[["Annotations"]]
  gene <- Gene
  a <- annot_df[annot_df[,2]==gene ,1]

  if(Normalized == TRUE) {

    count_data <- ExpCellNetObj[["Normalized Count Matrix"]]
    b <- mean(as.numeric(count_data[rownames(count_data) == a ,]))
    return(round(b , 1))
  }


  if( Normalized == FALSE) {

    count_data <- ExpCellNetObj[["Count Matrix"]]
    b <- mean(as.numeric(count_data[rownames(count_data) == a ,]))
    return(round(b , 1))
  }




}

########## Features Extraction ##########

#' Extract network-derived expression features
#'
#' Extracts a feature matrix from the normalized expression data by restricting
#' genes to those present in the constructed interaction network. The resulting
#' matrix is intended for downstream modeling and dimensionality reduction.
#'
#' @param ExpCellNetObj An ExpressionCellNet object containing a normalized count
#' matrix and a constructed network.
#'
#' @return A data.frame or matrix of network-derived expression features
#' (genes × samples), suitable for regression or multivariate analysis.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' features <- ExtractFeature(obj)
#' }

ExtractFeature <- function(ExpCellNetObj) {

  df <- ExpCellNetObj[["Network"]]
  annot_df <- ExpCellNetObj[["Annotations"]]
  Count_matrix <- ExpCellNetObj[["Normalized Count Matrix"]]


  all_genes <- c(df[,1] , df[,2])
  all_genes <- unique(all_genes)

  all_genes_ens <- annot_df[annot_df[,2] %in% all_genes,1]

  net_matrix <- Count_matrix[rownames(Count_matrix) %in% all_genes_ens,]

  for(i in 1:nrow(net_matrix)) {

    a <- rownames(net_matrix)[i]
    b <- annot_df[annot_df[,1] == a,2]
    rownames(net_matrix)[i] <- b

  }

  net_matrix <- t(net_matrix)

  Gen_df <- ExpCellNetObj[["Generations"]]


  net_matrix_sorted <- c()
  for(i in 1:nrow(Gen_df)) {
    a <- Gen_df[i,1]
    net_matrix_sorted <- cbind(net_matrix_sorted ,net_matrix[,colnames(net_matrix)== a] )
    colnames(net_matrix_sorted)[i] <- a
  }

  net_matrix_sorted <- as.data.frame(net_matrix_sorted)


  # remove Duplicates SCAMP3.1

  net_matrix_sorted <- t(net_matrix_sorted)
  net_matrix_sorted <- unique(net_matrix_sorted)
  net_matrix_sorted <- t(net_matrix_sorted)


  return(net_matrix_sorted)

}

########## Plot Correlation ##########

#' Plot correlation between two genes
#'
#' Plots correlation/regression between two genes (e.g., seed and distal gene).
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#' @param Gene1 Character; first gene.
#' @param Gene2 Character; second gene.
#'
#' @return Invisibly returns NULL (plot is produced).
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs
#'
#' @examples
#' \dontrun{
#' PlotCorrelation(obj, "PYGO2", "ATAD5")
#' }

PlotCorrelation <- function(ExpCellNetObj , Gene1 , Gene2) {
  all_feature <- ExtractFeature(ExpCellNetObj)

  all_feature <- all_feature[, colnames(all_feature) %in% c(Gene1, Gene2)]

  p1 <- ggplot2::ggplot(
    all_feature,
    ggplot2::aes(
      x = all_feature[, 1],
      y = all_feature[, 2]
    )
  ) +
    ggplot2::geom_point(color = "purple4") +
    ggplot2::labs(
      x = colnames(all_feature)[1],
      y = colnames(all_feature)[2],
      title = paste0(
        "Correlation Between ",
        colnames(all_feature)[1],
        " & ",
        colnames(all_feature)[2]
      )
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      se = FALSE,
      color = "tomato"
    )

  return(p1)


}

########## Regenerate Train Test Samples ##########

#' Regenerate train/test split
#'
#' Re-creates the internal random 80/20 train-test split without stratification.
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#'
#' @return Updated ExpressionCellNet object with regenerated split.
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- RegenerateTrainTest(obj)
#' }

RegenerateTrainTest <- function(ExpCellNetObj) {

  count_matrix <- ExpCellNetObj[["Count Matrix"]]

  train_number <- round(ncol(count_matrix) * 0.8 , 0)
  test_number <- ncol(count_matrix) - train_number

  train_samples <- sample(1:ncol(count_matrix) , train_number , replace = F)
  test_samples <- setdiff(1:ncol(count_matrix) , train_samples )

  train_df <- data.frame(col1 = colnames(count_matrix)[train_samples] , col2 = train_samples , col3 = rep("Train" , train_number))
  test_df <- data.frame(col1 = colnames(count_matrix)[test_samples] , col2 = test_samples , col3 = rep("Test" , test_number))

  Train_test_split <- rbind(train_df , test_df)

  colnames(Train_test_split) <- c("Sample_Name" , "Sample_Number" , "Group")

  Train_test_split <- as.data.frame(Train_test_split)

  if (length(ExpCellNetObj[["Train-Test Samples"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(Train_test_split))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Train-Test Samples"
    return(ExpCellNetObj)
  }

  if (length(ExpCellNetObj[["Train-Test Samples"]]) != 0){

    ExpCellNetObj[["Train-Test Samples"]] <- Train_test_split

    return(ExpCellNetObj)
  }




}


########## Load Train Test Split Data ##########
#' Load a predefined train/test split
#'
#' Replaces (or creates) the internal "Train-Test Samples" table in the object using a user-supplied split.
#' The input data.frame must include at least: Sample_Name, Sample_Number, Group (Train/Test).
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#' @param TrainTestSplitData A data.frame defining sample split into Train/Test.
#'
#' @return Updated ExpressionCellNet object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' split_df <- data.frame(
#'   Sample_Name = colnames(obj[["Count Matrix"]]),
#'   Sample_Number = seq_len(ncol(obj[["Count Matrix"]])),
#'   Group = c(rep("Train", 8), rep("Test", 2))
#' )
#' obj <- LoadTrainTestSplit(obj, split_df)
#' }

LoadTrainTestSplit <- function(ExpCellNetObj , TrainTestSplitData) {


  TrainTestSplitData <- as.data.frame(TrainTestSplitData)

  if (length(ExpCellNetObj[["Train-Test Samples"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(TrainTestSplitData))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Train-Test Samples"
    return(ExpCellNetObj)
  }

  if (length(ExpCellNetObj[["Train-Test Samples"]]) != 0){

    ExpCellNetObj[["Train-Test Samples"]] <- TrainTestSplitData

    return(ExpCellNetObj)
  }


}


########## Load Network Data ##########
#' Load network edge table
#'
#' Stores an externally provided edge list as the "Network" component of the object.
#' The table should include at least two columns corresponding to interacting gene symbols.
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#' @param NetworkData A data.frame representing network edges (two columns: gene1, gene2).
#'
#' @return Updated ExpressionCellNet object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' net <- data.frame(gene1=c("PYGO2","BRCA1"), gene2=c("UBAP2L","RAD51"))
#' obj <- LoadNetwork(obj, net)
#' }
#'
LoadNetwork <- function(ExpCellNetObj , NetworkData) {


  NetworkData <- as.data.frame(NetworkData)

  if (length(ExpCellNetObj[["Network"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(NetworkData))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Network"
    return(ExpCellNetObj)
  }

  if (length(ExpCellNetObj[["Network"]]) != 0){

    ExpCellNetObj[["Network"]] <- NetworkData

    return(ExpCellNetObj)
  }


}


########## Load Generations Data ##########
#' Load generations table
#'
#' Stores an externally provided gene-to-generation mapping as the "Generations" component.
#' The table should include columns Gene and Generation.
#'
#' @param ExpCellNetObj An ExpressionCellNet object.
#' @param GenerationsData A data.frame mapping genes to their generation labels.
#'
#' @return Updated ExpressionCellNet object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gen <- data.frame(Gene=c("PYGO2","UBAP2L"), Generation=c("Generation : 0","Generation : 1"))
#' obj <- LoadGenerations(obj, gen)
#' }
LoadGenerations <- function(ExpCellNetObj , GenerationsData) {


  GenerationsData <- as.data.frame(GenerationsData)

  if (length(ExpCellNetObj[["Generations"]]) == 0){

    ExpCellNetObj <- append(ExpCellNetObj , list(GenerationsData))
    names(ExpCellNetObj)[length(ExpCellNetObj)] <- "Generations"
    return(ExpCellNetObj)
  }

  if (length(ExpCellNetObj[["Generations"]]) != 0){

    ExpCellNetObj[["Generations"]] <- GenerationsData

    return(ExpCellNetObj)
  }


}






