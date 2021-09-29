library(shinydashboard)

# JS code to detect browser ----------------------------------------------------
js <- "
// execute the code after the shiny session has started
$(document).on('shiny:sessioninitialized', function(event) {
  // browser detection from https://stackoverflow.com/a/5918791/8099834
  navigator.sayswho= (function(){
    var ua= navigator.userAgent, tem, 
    M= ua.match(/(opera|chrome|safari|firefox|msie|trident(?=\\/))\\/?\\s*(\\d+)/i) || [];
    if(/trident/i.test(M[1])){
        tem=  /\\brv[ :]+(\\d+)/g.exec(ua) || [];
        return 'IE '+(tem[1] || '');
    }
    if(M[1]=== 'Chrome'){
        tem= ua.match(/\\b(OPR|Edge)\\/(\\d+)/);
        if(tem!= null) return tem.slice(1).join(' ').replace('OPR', 'Opera');
    }
    M= M[2]? [M[1], M[2]]: [navigator.appName, navigator.appVersion, '-?'];
    if((tem= ua.match(/version\\/(\\d+)/i))!= null) M.splice(1, 1, tem[1]);
    return M.join(' ');
  })(); 
  // pass browser info from JS to R
  Shiny.onInputChange('myBrowser', navigator.sayswho); 
});
"

tags$script(HTML(js))

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
    menuItem("String interaction tissue",
             tabName= "String-T",
             icon = icon("random", lib="glyphicon")),
    menuItem("String interaction blood",
             tabName= "String-B",
             icon = icon("tint", lib="glyphicon"))
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
            fluidRow(
              box(
                title = "Results table",
                width = 6,
                DTOutput("table_MA_T")
              ),
              box(title = "Forest plot",
                  plotlyOutput("forest_MA_T")
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
            fluidRow(
              box(
                title = "Results table",
                width = 6,
                DTOutput("table_MA_B")
              ),
              box(title = "Forest plot",
                  plotlyOutput("forest_MA_B")
              )
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
    tabItem(tabName = "String-T",
            p("This page shows the protein-protein interaction network computed
              by String database. The cutoff values are logFC=0.6, FDR=0.01 and
              a number of studies equal or greater than 15. To activate the 
              horizontal scrolling bar, please drag to the bottom the outter 
              vertical scrolling bar at the right of the page."),
            htmlOutput("stringT")
            ),
    tabItem(tabName = "String-B",
            p("This page shows the protein-protein interaction network computed
              by String database. The cutoff values are logFC=0.4, FDR=0.05 and
              a number of studies equal or greater than 2. To activate the 
              horizontal scrolling bar, please drag to the bottom the outter 
              vertical scrolling bar at the right of the page."),
            htmlOutput("stringB")
            )
  )
)

dashboardPage(header, sidebar, body)