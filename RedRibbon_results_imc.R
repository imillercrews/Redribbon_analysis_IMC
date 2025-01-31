#### Red ribbon functions 
# https://github.com/antpiron/RedRibbon

#### libraries ####
## install libraries
devtools::install_github("antpiron/RedRibbon")

# load libraries
library(RedRibbon)
library(tidyverse)

#### Create functions for Redribbon ####
### This melt.matrix function needed to be set up separetly for some reason
melt.matrix <- function (data, ..., na.rm = FALSE, value.name = "value")
{
  Var1 <- Var2 <- NULL
  
  dt <- as.data.table(data)
  colnames(dt) <- as.character(1:ncol(dt))
  dt[, rownames := 1:nrow(dt)]
  
  melted_dt <- data.table::melt(dt, id.vars = "rownames", na.rm = na.rm, value.name = value.name)
  colnames(melted_dt)  <- c("Var1", "Var2", "value")
  melted_dt[, Var1 := as.double(Var1)]
  melted_dt[, Var2 := as.double(Var2)]
  
  
  
  return(melted_dt)
}

### add value to set as max log scale for colors
#' Render the RRHO map.
#' 
#' You can choose to render the RedRibbon level map using \code{ggplot2}. 
#' 
#' @param self is the RedRibbon object
#' @param n is the number of coordinates to compute on the x and y axis (Default = sqrt(len))
#' @param labels is a list of labels of list a and b. Default is c("a", "b").
#' @param show.quadrants is a flag if set the quadrants lines are drawn
#' @param quadrants is the object returned by 'quadrants()' method
#' @param show.pval is a flag to show the P-values on the plot
#' @param repel.force is the value of the repel force for the p-value ploting (default = 150)
#' @param base_size is the size of the text fields (default = 20)
#' @param .log10 output log10 pval (default = FALSE)
#' @param ... The rest of the parameters
#' 
#' @return A \code{ggplot} object.
#' 
#' @method ggRedRibbon rrho
#' @export
ggRedRibbon.rrho.scale <- function (self, n = NULL, labels = c("a", "b"), show.quadrants=TRUE, quadrants=NULL, 
                                    show.pval=TRUE,
                                    repel.force=150, base_size=20, .log10=FALSE,
                                    new.max.log, # add value
                                    ...)
{
  len <- length(self$data$a)
  
  if ( is.null(n) )
    n <- min(max(sqrt(len), 500), len)
  
  n.i <- n
  n.j <- n
  
  rrho <- rrho_rectangle(1, 1, len, len, n.i, n.j, self$data$a, self$data$b,  mode=self$enrichment_mode, LOG=TRUE)
  if (.log10)
  {
    rrho <- rrho / log(10)
  }
  log.label <- if (.log10) "log10" else "log"
  
  # set top of log scale
  max.log <- new.max.log
  
  
  if (0 == max.log)
  {
    min.log  <- -0.001
    max.log <- 0.001
  } else
    min.log <- - max.log
  
  # remove negative p-value scale
  ticks <- c(min.log,0, max.log)
  
  len.colors <- length(self$ggplot_colours)
  half.len.colors <- len.colors %/% 2
  colors.values <- seq(0, len.colors) /  len.colors
  
  
  ## Suppress warning RRHO: no visible binding for global variable ‘gg’
  Var1 <- Var2 <- value <- i <- j <- pvalue <- NULL
  
  gg <-  ggplot2::ggplot(melt.matrix(rrho), ggplot2::aes(Var1,Var2, fill=value)) +
    ggplot2::geom_raster() +
    ## ggplot2::scale_fill_gradientn(colours=self$ggplot_colours, name="-log p.val") +
    ggplot2::scale_fill_gradientn(colors = self$ggplot_colours,
                                  breaks = ticks,
                                  labels = format(ticks),
                                  limits=ticks[c(1,3)],
                                  ##limits=b[c(1,length(colors))],
                                  name=paste0("-", log.label, " p.val"),
                                  values=colors.values) +
    ggplot2::xlab(labels[1]) + ggplot2::ylab(labels[2]) +
    ## scale_x_continuous(labels = label_percent(accuracy = 1, scale = 100/n.i)) +
    ## scale_y_continuous(labels = label_percent(accuracy = 1, scale = 100/n.j) ) +
    ggplot2::scale_x_continuous(breaks = c(0 + n * 0.1, n - n * 0.1), labels = c("down", "up"), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = c(0 + n * 0.1, n - n * 0.1), labels = c("down", "up"), expand = c(0, 0)) +
    ## ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_text(size=base_size,face="bold"),
                   legend.title = ggplot2::element_text(size = base_size * 7 / 10),
                   legend.text = ggplot2::element_text(size = base_size * 1 / 2),
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=base_size* 7 / 10, face="bold"),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size=base_size * 7 / 10, face="bold", angle=90),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.ticks.length = ggplot2::unit(0, "pt"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(), 
                   panel.spacing = ggplot2::unit(0, "cm"),
                   plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"))
  
  
  ## find the middle of the plots
  if (show.quadrants || show.pval)
  {
    a_ltzero <- sum(self$data$a < 0)
    x.ind <- a_ltzero
    if ( x.ind == 0 )
      x.ind = len / 2
    
    b_ltzero <- sum(self$data$b < 0)
    y.ind <- b_ltzero
    if ( y.ind == 0 )
      y.ind = len / 2
    
    ## plot dotted quadrant lines
    if (show.quadrants)
    {
      gg  <- gg +
        ggplot2::geom_vline(ggplot2::aes(xintercept = x.ind * n.i / len), 
                            linetype = "dotted", colour = "gray10",size = 1) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = y.ind * n.j / len), 
                            linetype = "dotted", colour = "gray10",size = 1)
    }
    
    ## plot pvalue
    if (show.pval)
    {
      if (! is.null(quadrants) )
      {
        pval_size  <- as.integer(base_size * 1/5)
        quadrants_df <- as.data.frame(
          do.call(rbind, lapply(quadrants,
                                function (quadrant)
                                {
                                  if ( quadrant$pvalue > 0.05 || (! is.null(quadrant$padj) && quadrant$padj > 0.05) )
                                    return(NULL)
                                  
                                  pvalue <- if (.log10) quadrant$log_pvalue / log(10) else quadrant$log_pvalue
                                  pvalue.formatted <-  formatC(pvalue, format = "f", digits = 1)
                                  
                                  if (! is.null(quadrant$padj) )
                                  {
                                    padj <- if (.log10) quadrant$log_padj / log(10) else quadrant$log_padj
                                    pvalue.formatted <-  paste(pvalue.formatted,
                                                               "(padj =", formatC(padj, format = "f", digits = 1), ")")
                                  }
                                  data.frame(i=quadrant$i, j=quadrant$j,
                                             pvalue=pvalue.formatted, value=pvalue)
                                })))
        
        if ( nrow(quadrants_df) > 0 )
          gg <- gg +
          ggrepel::geom_text_repel(data=quadrants_df,
                                   ggplot2::aes(x=i * n.i / len, y=j * n.j / len,
                                                label=pvalue,
                                                colour = "gray"),
                                   hjust=1, vjust=1, colour = "black",
                                   force = repel.force, show.legend = FALSE, size = pval_size)
        
      }
    }
  }
  
  return(gg)
}

#### create function to output redribbon graph and save quadrant data
### just input both gene lists
#' @param celltype is the title of the analysis
#' @param dataset.a is the DEG gene list with one column the gene names and the second column the logFC or p.value results with positive and negative indicating directionality 
#' @param dataset.b is the same as dataset.a
#' @param dataset.a.type is the name given to the 'a' gene list analysis
#' @param dataset.b.type is the same as @param dataset.a.type
#' @param a.variable is the variable name with the values from dataset.a
#' @param b.variable is the same as a.variable
#' @param new.max.log is the value to set the max and min for the log scale in the heatmap graph. Can be left 'NULL' if scale should be determined by the redribbon function. The new.max.log needs to be at least the value of the log scale from the redribbon function. Useful if you want to set the scale to the max value across mutlitple redribbon objects
#' @param file.name is the folder file path where the figure and quadrant data should be saved
RedRibbon.all <- function(celltype = "Analysis.title",
                          dataset.a = gene.list.a,
                          dataset.b = gene.list.b,
                          dataset.a.type = "Name.a",
                          dataset.b.type = "Name.b",
                          a.variable = "Value.name.a",
                          b.variable = "Value.name.b",
                          new.max.log = NULL,
                          file.name = './file.name/')
{
  # create data frame for red ribbon
  # needs to have an id (gene) col and one called 'a' and one called 'b'
  df = dataset.a %>% 
    rename('a' =  a.variable) %>% 
    full_join(dataset.b %>% 
                rename('b' = b.variable))
  
  ## Create RedRibbon object
  rr <- RedRibbon(df, 
                  enrichment_mode="hyper-two-tailed")
  
  ## Run the overlap using evolutionnary algorithm,
  ## computing permutation adjusted p-value for the four quadrants
  quad <- quadrants(rr, 
                    algorithm="ea",
                    permutation=TRUE, 
                    whole=FALSE)
  ### compare RRHO2 to Redribbon
  # create list of RRHO outcomes
  RR.list <- data.frame(Gene = df[quad$upup$positions,]$gene,
                        RRquadrant = 'upup') %>% 
    rbind(data.frame(Gene = df[quad$downdown$positions,]$gene,
                     RRquadrant = 'downdown')) %>% 
    rbind(data.frame(Gene = df[quad$updown$positions,]$gene,
                     RRquadrant = 'updown')) %>% 
    rbind(data.frame(Gene = df[quad$downup$positions,]$gene,
                     RRquadrant = 'downup')) %>% 
    mutate(Sample = celltype)
  
  # save file
  write_csv(RR.list,
            file = paste0(file.name,
                          celltype,
                          ' quadrant genes.csv'))
  ## Plots the results
  ggRedRibbon(rr,
              quadrants=quad) + 
    coord_fixed(ratio = 1, 
                clip = "off") +
    xlab(dataset.a) +
    ylab(dataset.b) +
    ggtitle(celltype)
  
  if (is.null(new.max.log))
  {
    ggRedRibbon(rr, 
                quadrants=quad) + 
      coord_fixed(ratio = 1, 
                  clip = "off") +
      xlab(dataset.a.type) +
      ylab(dataset.b.type) +
      ggtitle(celltype)
    ggsave(paste0(file.name,
                  celltype,
                  ' RedRibbon.png'),
           height = 10,
           width = 10)
  } else
  {
  # scaled
  ggRedRibbon.rrho.scale(rr, 
                         quadrants=quad,
                         new.max.log = new.max.log) + 
    coord_fixed(ratio = 1, 
                clip = "off") +
    xlab(dataset.a.type) +
    ylab(dataset.b.type) +
    ggtitle(celltype)
  ggsave(paste0(file.name,
                celltype,
                ' RedRibbon_scaled.png'),
         height = 10,
         width = 10)
  }

  
}

#### Run Redribbon ####
### parameter description for new function

#' @param celltype is the title of the analysis
#' @param dataset.a is the DEG gene list with one column the gene names and the second column the logFC or p.value results with positive and negative indicating directionality 
#' @param dataset.b is the same as dataset.a
#' @param dataset.a.type is the name given to the 'a' gene list analysis
#' @param dataset.b.type is the same as @param dataset.a.type
#' @param a.variable is the variable name with the values from dataset.a
#' @param b.variable is the same as a.variable
#' @param new.max.log is the value to set the max and min for the log scale in the heatmap graph. Can be left 'NULL' if scale should be determined by the redribbon function. The new.max.log needs to be at least the value of the log scale from the redribbon function. Useful if you want to set the scale to the max value across mutlitple redribbon objects
#' @param file.name is the folder file path where the figure and quadrant data should be saved


### Run function

RedRibbon.all(celltype = "Analysis.title",
                          dataset.a = gene.list.a,
                          dataset.b = gene.list.b,
                          dataset.a.type = "Name.a",
                          dataset.b.type = "Name.b",
                          a.variable = "Value.name.a",
                          b.variable = "Value.name.b",
                          new.max.log = NULL,
                          file.name = './file.name/')
