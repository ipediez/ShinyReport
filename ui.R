library(shinydashboard)

header <- dashboardHeader(title = "IVO-UBB: PDAC")
#header <- dashboardHeader(title = span(img(src="UBB.png", width = 190)))

sidebar <- dashboardSidebar(
  img(src="CIPF.png", width = "100%"),
  img(src="UBB.png", width = "100%"),
  sidebarMenu(
    menuItem("Differential Expression - Tissue", 
             tabName = "DE-tissue",
             icon = icon("equalizer", lib="glyphicon")),
    menuItem("Differential Expression - Blood", 
             tabName = "DE-blood",
             icon = icon("tint", lib="glyphicon")),
    menuItem("DE Meta-analysis - Tissue", 
             tabName = "MA-tissue",
             icon = icon("stats", lib="glyphicon")),
    menuItem("DE Meta-analysis - Blood", 
             tabName = "MA-blood",
             icon = icon("tint", lib="glyphicon")),
    menuItem("DE Meta-analysis \n Tissue vs Blood",
             tabName = "Tissue-Blood",
             icon = icon("magnet", lib="glyphicon")),
    menuItem(a("String interaction tissue", href="https://version-11-5.string-db.org/cgi/network?networkId=bhxSHXGN5uzV", target="_blank"),
             tabName= "String-T"),
    menuItem(a("String interaction blood", href="https://version-11-5.string-db.org/cgi/network?networkId=bhxSHXGN5uzV", target="_blank"),
             tabName= "String-B"),
    menuItem("Gene Ontology \n Biological Process",
             tabName = "GO-BP",
             icon = icon("stats", lib="glyphicon"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "DE-tissue",
            p("El análisis de expresión diferencial se ha realizado utilizando
              el paquete de R 'limma' para los datos de microarray, y el paquete
              de R 'DESeq2' para los datos de RNAseq. Se han tenido en cuenta
              los efectos batch y las muestras emparejadas."),
            fluidRow(
              tabBox(
                title = "Summary",
                id = "De-tissue-tabset1", height = "250px",
                tabPanel("Stats",
                         textOutput("total_genes_DE_T"),
                         textOutput("pval_genes_DE_T"),
                         textOutput("up_genes_DE_T"),
                         textOutput("down_genes_DE_T")),
                tabPanel("Vulcano plot",
                         plotlyOutput("vulcano_DE_T",
                                      height=350))
              ),
              box(title = "Inputs",
                  selectInput("dataset_tissue",
                              label = h3("Select a dataset"),
                              choices = names(expression_tissue), 
                              selected = names(expression_tissue)[1]),
                  numericInput("logfc_DE_T", h3("Select logFC threshold"),
                               value=0.6),
                  numericInput("FDR_DE_T", h3("Select FDR threshold"),
                               value=0.05))
            ),
            downloadButton('download_DE_T', "Download the data"),
            fluidRow(
              box(
                title = "Results table",
                width = 12,
                DTOutput("table_DE_T")
              )
            )),
    tabItem(tabName = "DE-blood",
            p("El análisis de expresión diferencial se ha realizado utilizando
              el paquete de R 'limma' para los datos de microarray. Se han
              tenido en cuenta los efectos batch y las muestras emparejadas."),
            fluidRow(
              tabBox(
                title = "Summary",
                id = "De-blood-tabset1", height = "250px",
                tabPanel("Stats",
                         textOutput("total_genes_DE_B"),
                         textOutput("pval_genes_DE_B"),
                         textOutput("up_genes_DE_B"),
                         textOutput("down_genes_DE_B")),
                tabPanel("Vulcano plot",
                         plotlyOutput("vulcano_DE_B",
                                      height=350))
              ),
              box(title = "Inputs",
                  selectInput("dataset_blood",
                              label = h3("Select a dataset"),
                              choices = names(expression_blood), 
                              selected = names(expression_blood)[1]),
                  numericInput("logfc_DE_B", h3("Select logFC threshold"),
                               value=0.6),
                  numericInput("FDR_DE_B", h3("Select FDR threshold"),
                               value=0.05))
            ),
            downloadButton('download_DE_B', "Download the data"),
            fluidRow(
              box(
                title = "Results table",
                width = 12,
                DTOutput("table_DE_B")
              )
            )),
    tabItem(tabName = "MA-tissue",
            p("Meta-analysis has been performed following Dersimonian-Laird methond
              with the R package 'metafor'."),
            fluidRow(
              tabBox(
                title = "Summary",
                id = "Ma-tissue-tabset1", height = "250px",
                tabPanel("Stats",
                         textOutput("total_genes_MA_T"),
                         textOutput("pval_genes_MA_T"),
                         textOutput("up_genes_MA_T"),
                         textOutput("down_genes_MA_T")),
                tabPanel("Vulcano plot",
                         plotlyOutput("vulcano_MA_T",
                                      height=450))
              ),
              box(title = "Inputs",
                  numericInput("logfc_MA_T", h3("Select logFC threshold:"),
                               value=0.6),
                  numericInput("FDR_MA_T", h3("Select FDR threshold:"),
                               value=0.05),
                  numericInput("studies_MA_T", h3("Select minimum number of studies:"),
                               value=2),
                  textInput("Gene_MA_T", h3("Forest plot:"),
                            "Type an ENTREZ ID"))
            ),
            downloadButton('download_MA_T', "Download the data"),
            fluidRow(
              box(
                title = "Results table",
                width = 6,
                DTOutput("table_MA_T")
              ),
              conditionalPanel(condition ="input.Gene_MA_T != 'Type an ENTREZ ID'",
              box(title = "Forest plot",
                  plotlyOutput("forest_MA_T"))
              )
            )),
    tabItem(tabName = "MA-blood",
            p("Meta-analysis has been performed following Dersimonian-Laird methond
              with the R package 'metafor'."),
            fluidRow(
              tabBox(
                title = "Summary",
                id = "Ma-blood-tabset1", height = "250px",
                tabPanel("Stats",
                         textOutput("total_genes_MA_B"),
                         textOutput("pval_genes_MA_B"),
                         textOutput("up_genes_MA_B"),
                         textOutput("down_genes_MA_B")),
                tabPanel("Vulcano plot",
                         plotlyOutput("vulcano_MA_B",
                                      height=450))
              ),
              box(title = "Inputs",
                  numericInput("logfc_MA_B", h3("Select logFC threshold:"),
                               value=0.6),
                  numericInput("FDR_MA_B", h3("Select FDR threshold:"),
                               value=0.05),
                  numericInput("studies_MA_B", h3("Select minimum number of studies:"),
                               value=2),
                  textInput("Gene_MA_B", h3("Forest plot:"),
                            "Type an ENTREZ ID"))
            ),
            downloadButton('download_MA_B', "Download the data"),
            fluidRow(
              box(
                title = "Results table",
                width = 6,
                DTOutput("table_MA_B")
              ),
              conditionalPanel(condition ="input.Gene_MA_B != 'Type an ENTREZ ID'",
              box(title = "Forest plot",
                  plotlyOutput("forest_MA_B")
              ))
            )),
    tabItem(tabName = "Tissue-Blood",
            p("This page shows the intersection of genes between both meta-analyses.
              The values for logFC, FDR and number of studies used to decide
              wherther a gene is up/down regulated or not are taken from the 
              previous two pages."),
            fluidRow(
              box(
                title = "Both upregulated",
                DTOutput("compare_UP")
              ), 
              box(
                title = "Both downregulated",
                DTOutput("compare_DOWN")
              )
            ),
            fluidRow(
              box(
                title = "Tissue up, blood down",
                DTOutput("compare_UD")
              ),
              box(
                title = "Tissue down, blood up",
                DTOutput("compare_DU")
              )
            )
      
    ),
    tabItem(tabName = "GO-BP",
            p("This page shows the results of functional enrichment with MDGSA R
            package. The followin filters have been applied: adjusted p value
              lower than 0.05, functions are mapped to at least 11 genes and
              less than 200 genes. This last filter avoids too specific ontologies
              and too wide ones."),
            fluidRow(
              box(
                width = 12,
                column(4,
                       # Copy the line below to make a slider range 
                       sliderInput("slider_bp", label = h3("Genes Range"), min = 10, 
                                   max = 200, value = c(10, 200)))
              )
            ),
            fluidRow(
              box(
                title = "Gene Ontology: Biological Process",
                width = 12,
                DTOutput("enrichment_BP")
              ))
            
    )
  )
)

dashboardPage(header, sidebar, body)