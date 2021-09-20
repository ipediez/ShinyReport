library(shiny)
library(plotly)
library(tidyverse)
library(kableExtra)
library(DT)

load("DEdata.Rdata")

# Define Vulcano function ----

vulcano_plotly <- function(input){
  data <- expression[[input$dataset]]
  
  data$diffexpressed <- "NO"
  data$diffexpressed[data$logFC > input$logfc & data$adj.P.Val < input$FDR] <- "UP"
  data$diffexpressed[data$logFC < -input$logfc & data$adj.P.Val < input$FDR] <- "DOWN"
  
  
  data$ENTREZID <- rownames(data)
  
  fig <- plot_ly(data, x=~logFC, y=~log10(adj.P.Val),
                 type="scatter",
                 mode="markers",
                 color=~diffexpressed, colors= c("#2166ac", "black", "#b2182b"),
                 text= ~paste("ENTREZ ID: ", ENTREZID, "<br>Gene name: ", geneNames)
  ) %>%
    layout(
      yaxis=list(autorange="reversed")
    )
  return(fig)
}

# Define table function ----
my_table <- function(input){
  data <- expression[[input$dataset]]
  
  data$diffexpressed <- "NO"
  data$diffexpressed[data$logFC > input$logfc & data$adj.P.Val < input$FDR] <- "UP"
  data$diffexpressed[data$logFC < -input$logfc & data$adj.P.Val < input$FDR] <- "DOWN"
  
  dt_link <- sapply(rownames(data), function(x) paste0("https://www.ncbi.nlm.nih.gov/gene/", x))
  
  my_tab <- data %>%
    mutate(geneNames= cell_spec(geneNames, "html", link = dt_link)) %>%
    filter(diffexpressed %in% c("UP", "DOWN")) %>%
    select(logFC, se.coef, P.Value, adj.P.Val, geneNames) %>%
    datatable(escape=FALSE, filter = "top", options=list(pageLength=10, autoWidth=TRUE))
  
  return(my_tab)
}

# Define summary function ----
sum_values <- function(input){
  data <- expression[[input$dataset]]
  
  data$diffexpressed <- "NO"
  data$diffexpressed[data$logFC > input$logfc & data$adj.P.Val < input$FDR] <- "UP"
  data$diffexpressed[data$logFC < -input$logfc & data$adj.P.Val < input$FDR] <- "DOWN"
  
  total_genes <- dim(data)[1]
  pval_genes <- length(data$diffexpressed[data$adj.P.Val < input$FDR])
  up_genes <- dim(data[data$diffexpressed == "UP",])[1]
  down_genes <- dim(data[data$diffexpressed == "DOWN",])[1]
  
  return(c(total_genes, pval_genes, up_genes, down_genes))
}


# Define UI ----

ui <- fluidPage(
  titlePanel("Differential Expression"),
  
  sidebarLayout(
    sidebarPanel(
      htmlOutput("UBB"),
      selectInput("dataset",
                  label = h3("Select a dataset"),
                  choices = names(expression), 
                  selected = names(expression)[1]),
      numericInput("logfc", h3("Select logFC threshold"),
                   value=0.6),
      numericInput("FDR", h3("Select FDR threshold"),
                   value=0.05)
    ),
    mainPanel(p("El análisis de expresión diferencial se ha realizado utilizando el paquete de R limma para los estudios de microarray y DESeq2 para los estudios de RNA-seq"),
              h1("Resumen"),
              textOutput("total_genes"),
              textOutput("pval_genes"),
              textOutput("up_genes"),
              textOutput("down_genes"),
              h1("Vulcano Plot"),
              plotlyOutput("vulcano"),
              h1("Tabla de resultados"),
              DTOutput("table"))
  ),
  

)

# Define server logic ----

server <- function(input, output) {
  output$vulcano <- renderPlotly(
    vulcano_plotly(input)
  )
  
  output$table <- renderDataTable(
    my_table(input)
  )
  
  output$total_genes <- renderText(
    paste("Numero total de genes analizados:", sum_values(input)[1])
  )
  
  output$pval_genes <- renderText(
    paste0("Número de genes significativos (FDR < ", input$FDR, "): ", sum_values(input)[2])
  )
  
  output$up_genes <- renderText(
    paste("Número de genes UP regulados (logFC > ", input$logfc, "): ", sum_values(input)[3])
  )
  
  output$down_genes <- renderText(
    paste("Numero de genes DOWN regulados (logFC < -", input$logfc, "): ", sum_values(input)[4])
  )
  
  output$UBB <- renderText({c('<img src="',"https://bioinfo.cipf.es/ubb/wp-content/uploads/2018/03/logos.png",'">')})
}

# Run the app ----

shinyApp(ui = ui, server = server)