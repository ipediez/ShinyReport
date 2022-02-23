library(shiny)
library(plotly)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(kableExtra)
library(DT)

# LOAD DATA --------------------------------------------------------------------
expression_tissue <- readRDS("./Data/DE_tissue.RDS")
expression_blood <- readRDS("./Data/DE_blood.RDS")
expression <- c(expression_tissue, expression_blood)

table_MA_T <- read.csv2("./Data/all.results.tissue.txt", sep="\t", stringsAsFactors = F)
results_MA_T <- readRDS("./Data/tissue_metaanalysis.RDS")
table_MA_B <- read.csv2("./Data/all.results.blood.txt", sep="\t", stringsAsFactors = F)
results_MA_B <- readRDS("./Data/blood_metaanalysis.RDS")

BP_enrichment <- read.csv("../Data/Functional_enrichment/mdgsa_tissue/sig_metaanalysis_bp.tsv", sep = "\t")

# COMMON FUNCTIONS -------------------------------------------------------------
filter_data <- function(data, logFC, FDR, studies=FALSE){
  if (studies == FALSE) {
    data$diffexpressed <- "NO"
    data$diffexpressed[data$logFC > logFC & data$adj.P.Val < FDR] <- "UP"
    data$diffexpressed[data$logFC < -logFC & data$adj.P.Val < FDR] <- "DOWN"
    
    data <- data %>% rename(FDR = adj.P.Val)
  } else {
    rownames(data) <- data$ID
    data$summary_logFC <- as.numeric(data$summary_logFC)
    data$p.adjust <- as.numeric(data$p.adjust)
    
    data$diffexpressed <- "NO"
    data$diffexpressed[data$summary_logFC > logFC &
                         data$p.adjust < FDR &
                         data$studies >= studies] <- "UP"
    data$diffexpressed[data$summary_logFC < -logFC &
                         data$p.adjust < FDR &
                         data$studies >= studies] <- "DOWN"
    
    data <- data %>% rename(FDR = p.adjust, logFC = summary_logFC, geneNames = name)
    
  }
  return(data)
}

## Define Vulcano plot function------------------------------------------------
vulcano <- function(data, logFC, FDR, studies=FALSE){
  
  data <- filter_data(data, logFC, FDR, studies)
  
  col <- c("#2166ac", "black", "#b2182b")
  names(col) <- c("DOWN", "NO", "UP")
  data$ENTREZID <- rownames(data)
  
  fig <- plot_ly(data, x=~logFC, y=~log10(FDR),
                 type="scatter",
                 mode="markers",
                 color=~diffexpressed, colors= col,
                 opacity = "0.6",
                 text= ~paste("ENTREZ ID: ", ENTREZID, "<br>Gene name: ", geneNames)
  ) %>%
    layout(
      yaxis=list(autorange="reversed")
    )
  return(fig)
}

## Define summary function --------------------------------------------------
sum_values <- function(data, logFC, FDR, studies=FALSE){
  
  data <- filter_data(data, logFC, FDR, studies)
  
  total_genes <- dim(data)[1]
  pval_genes <- dim(data[data$diffexpressed %in% c("UP","DOWN"),])[1]
  up_genes <- dim(data[data$diffexpressed == "UP",])[1]
  down_genes <- dim(data[data$diffexpressed == "DOWN",])[1]
  
  return(c(total_genes, pval_genes, up_genes, down_genes))
}

## Define table function  ------------------------------------------------------
my_table <- function(data, logFC, FDR, studies=FALSE){
  if (studies == FALSE) {
    data$diffexpressed <- "NO"
    data$diffexpressed[data$logFC > logFC & data$adj.P.Val < FDR] <- "UP"
    data$diffexpressed[data$logFC < -logFC & data$adj.P.Val < FDR] <- "DOWN"
    
    dt_link <- sapply(rownames(data), function(x) paste0("https://www.ncbi.nlm.nih.gov/gene/", x))
    
    my_tab <- data %>%
      mutate(geneNames= cell_spec(geneNames, "html", link = dt_link)) %>%
      filter(diffexpressed %in% c("UP", "DOWN")) %>%
      select(logFC, se.coef, P.Value, adj.P.Val, geneNames) %>%
      datatable(escape=FALSE, filter = "top", options=list(pageLength=10, autoWidth=TRUE))
  } else{
    rownames(data) <- data$ID
    data$summary_logFC <- as.numeric(data$summary_logFC)
    data$p.adjust <- as.numeric(data$p.adjust)
    
    data$diffexpressed <- "NO"
    data$diffexpressed[data$summary_logFC > logFC &
                         data$p.adjust < FDR &
                         data$studies >= studies] <- "UP"
    data$diffexpressed[data$summary_logFC < -logFC &
                         data$p.adjust < FDR &
                         data$studies >= studies] <- "DOWN"
    
    dt_link <- sapply(data$ID, function(x) paste0("https://www.ncbi.nlm.nih.gov/gene/", x))
    
    my_tab <- data %>%
      mutate(name= cell_spec(name, "html", link = dt_link)) %>%
      filter(diffexpressed %in% c("UP", "DOWN")) %>%
      select(ID, summary_logFC, p.adjust, name, studies) %>%
      rename(logFC=summary_logFC, FDR=p.adjust) %>%
      datatable(escape=FALSE,
                filter = "top",
                options=list(pageLength=6,
                             autoWidth=TRUE),
                rownames = FALSE)
  }
  
  return(my_tab)
}


# META-ANALYSIS EXCLUSIVE FUNCTIONS --------------------------------------------
## Define forest plot ----------------------------------------------------------
forest_plot <- function(gene, data, meta_analysis){
  if (gene == "Type an ENTREZ ID") {
    return("No gene")
  } else {
    
    res <- meta_analysis[[as.character(gene)]]
    rownames(data) <- data$ID
    
    # Prepare data frame ---------------------------------------------------------
    overall_model <- c("Overall Model", res$b, res$ci.lb, res$ci.ub, "overall")
    
    # Compute confidence interval for the individual values of each study
    ci.lb <- res$yi - qnorm(res$level/2, lower.tail=FALSE) * sqrt(as.numeric(res$vi))
    ci.ub <- res$yi + qnorm(res$level/2, lower.tail=FALSE) * sqrt(as.numeric(res$vi))
    
    individual_values <- data.frame(names(res$yi),
                                    res$yi,
                                    ci.lb,
                                    ci.ub,
                                    rep("individual", length(res$yi)),
                                    stringsAsFactors = F)
    
    plot_data <- rbind(individual_values, overall_model)
    
    plot_data <- cbind(c(1:dim(plot_data)[1]), plot_data)
    colnames(plot_data) <- c("Index", "Study", "logFC", "ci_l", "ci_u", "highlight")
    plot_data$logFC <- as.numeric(plot_data$logFC)
    plot_data$ci_l <- as.numeric(plot_data$ci_l)
    plot_data$ci_u <- as.numeric(plot_data$ci_u)
    mycolours <- c("overall" = "red", "individual" = "black")
    
    labels <- ifelse(plot_data$logFC > 0,
                     paste0(" ",
                            round(plot_data$logFC, 2),
                            " [",
                            round(plot_data$ci_l, 2),
                            ", ",
                            round(plot_data$ci_u, 2),
                            "]"),
                     paste0(round(plot_data$logFC, 2),
                            " [",
                            round(plot_data$ci_l, 2),
                            ", ",
                            round(plot_data$ci_u, 2),
                            "]"))
    
    plot_data$labels <- labels
    
    # Plot -----------------------------------------------------------------------
    
    # )
    p <- ggplot(data=plot_data, aes(y=Index, x=logFC, xmin=ci_l, xmax=ci_u, label = labels))
    # adds the effect sizes to the plot
    p <- p + geom_point(aes(colour=highlight))
    p <- p + scale_color_manual("Status", values=mycolours)
    # adds the CIs
    p <- p + geom_errorbarh(height=.1, aes(colour=highlight))
    # add the labels
    p <- p + geom_text(mapping = aes(y=Index, x=3), hjust=-0.1, size=3)
    p <- p + scale_x_continuous(limits=c(-3,3), breaks = c(-3:3), name="logFC")
    p <- p + scale_y_continuous(name = "", breaks=1:dim(plot_data)[1], labels = plot_data$Study, trans="reverse")
    p <- p + geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)
    p <- p + ggpubr::theme_pubr()
    p <- p + theme(plot.margin = unit(c(1, 5, 1, 1), "lines"), legend.position = "none") +
      coord_cartesian(clip="off")
    title <- paste(as.character(gene), data[as.character(gene), "name"], sep = "\n")
    p <- p + ggtitle(title)
    
    fig <- ggplotly(p)
    
    return(fig)
  }
  
}

## Compare meta-analysis tables ------------------------------------------------

compare_table <- function(first_table, second_table, input, first_sign="UP",
                          second_sign="UP", names=c("tissue", "blood")){
  #Filter data
  tab1 <- filter_data(data = first_table,
                      logFC = input$logfc_MA_T,
                      FDR = input$FDR_MA_T,
                      studies = input$studies_MA_T)
  tab2 <- filter_data(data = second_table,
                      logFC = input$logfc_MA_B,
                      FDR = input$FDR_MA_B,
                      studies = input$studies_MA_B)
  
  newnames <- c(paste("logFC", names[1]), paste("logFC", names[2]))
  
  x <- left_join(tab1, tab2, by="ID") %>%
    filter(diffexpressed.x == first_sign, diffexpressed.y == second_sign) %>%
    mutate(name = cell_spec(geneNames.x, "html", link =paste0("https://www.ncbi.nlm.nih.gov/gene/", ID))) %>%
    select(ID, name, logFC.x, logFC.y) %>%
    rename_at(vars(c("logFC.x", "logFC.y")), ~ newnames) %>%
    datatable(escape=FALSE,
              filter = "top",
              options=list(pageLength= 10,
                           autoWidth=TRUE),
              rownames = FALSE)
    
  return(x)
}

shinyAppDir("./")

# Functional enrichment --------------------------------------------------------
## Define GO enrichment table
bp_table <- function(data){
  my_tab <- data %>%
    datatable(escape=FALSE,
              filter = "top",
              options=list(pageLength=6,
                           autoWidth=TRUE),
              rownames = TRUE)
  
  return(my_tab)
}