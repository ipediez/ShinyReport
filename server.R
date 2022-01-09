server <- function(input, output) {
  
  # Differential Expression - Tissue -------------------------------------------
  output$vulcano_DE_T <- renderPlotly(
    vulcano(expression[[input$dataset_tissue]], input$logfc_DE_T, input$FDR_DE_T)
  )
  
  output$table_DE_T <- renderDataTable(
    my_table(expression[[input$dataset_tissue]], input$logfc_DE_T, input$FDR_DE_T)
  )
  output$download_DE_T <- downloadHandler(
    filename = function(){paste0("DEA_", input$dataset_tissue, ".csv")},
    content = function(fname){
      write.csv(my_table(expression[[input$dataset_tissue]], input$logfc_DE_T, input$FDR_DE_T), fname)
    }
  )

  output$total_genes_DE_T <- renderText(
    paste("Numero total de genes analizados:",
          sum_values(expression[[input$dataset_tissue]], input$logfc_DE_T, input$FDR_DE_T)[1])
  )

  output$pval_genes_DE_T <- renderText(
    paste0("Número total de genes significativos (FDR < ", input$FDR_DE_T, "): ",
           sum_values(expression[[input$dataset_tissue]], input$logfc_DE_T, input$FDR_DE_T)[2])
  )

  output$up_genes_DE_T <- renderText(
    paste("Número de genes UP regulados (logFC > ", input$logfc_DE_T, "): ",
          sum_values(expression[[input$dataset_tissue]], input$logfc_DE_T, input$FDR_DE_T)[3])
  )

  output$down_genes_DE_T <- renderText(
    paste("Numero de genes DOWN regulados (logFC < -", input$logfc_DE_T, "): ",
          sum_values(expression[[input$dataset_tissue]], input$logfc_DE_T, input$FDR_DE_T)[4])
  )
  
  # Differential Expression - Blood -------------------------------------------
  output$vulcano_DE_B <- renderPlotly(
    vulcano(expression[[input$dataset_blood]], input$logfc_DE_B, input$FDR_DE_B)
  )
  
  output$table_DE_B <- renderDataTable(
    my_table(expression[[input$dataset_blood]], input$logfc_DE_B, input$FDR_DE_B)
  )
  
  output$download_DE_B <- downloadHandler(
    filename = function(){paste0("DEA_", input$dataset_tissue, ".csv")},
    content = function(fname){
      write.csv(my_table(expression[[input$dataset_blood]], input$logfc_DE_B, input$FDR_DE_B), fname)
    }
  )
  
  output$total_genes_DE_B <- renderText(
    paste("Numero total de genes analizados:",
          sum_values(expression[[input$dataset_blood]], input$logfc_DE_B, input$FDR_DE_B)[1])
  )
  
  output$pval_genes_DE_B <- renderText(
    paste0("Número total de genes significativos (FDR < ", input$FDR_DE_B, "): ",
           sum_values(expression[[input$dataset_blood]], input$logfc_DE_B, input$FDR_DE_B)[2])
  )
  
  output$up_genes_DE_B <- renderText(
    paste("Número de genes UP regulados (logFC > ", input$logfc_DE_B, "): ",
          sum_values(expression[[input$dataset_blood]], input$logfc_DE_B, input$FDR_DE_B)[3])
  )
  
  output$down_genes_DE_B <- renderText(
    paste("Numero de genes DOWN regulados (logFC < -", input$logfc_DE_B, "): ",
          sum_values(expression[[input$dataset_blood]], input$logfc_DE_B, input$FDR_DE_B)[4])
  )
  
  # Meta-analysis - Tissue -----------------------------------------------------
  output$vulcano_MA_T <- renderPlotly(
    vulcano(table_MA_T, input$logfc_MA_T, input$FDR_MA_T, input$studies_MA_T)
  )
  
  output$forest_MA_T <- renderPlotly(
    forest_plot(input$Gene_MA_T, data=table_MA_T, meta_analysis=results_MA_T)
  )
  
  output$table_MA_T <- renderDataTable(
    my_table(table_MA_T, input$logfc_MA_T, input$FDR_MA_T, input$studies_MA_T)
  )
  
  output$download_MA_T <- downloadHandler(
    filename = function(){"metaanalysis.csv"},
    content = function(fname){
      write.csv(my_table(table_MA_T, input$logfc_MA_T, input$FDR_MA_T, input$studies_MA_T), fname)
    }
  )
  
  output$total_genes_MA_T <- renderText(
    paste("Numero total de genes analizados:",
          sum_values(table_MA_T, input$logfc_MA_T, input$FDR_MA_T, input$studies_MA_T)[1])
  )
  
  output$pval_genes_MA_T <- renderText(
    paste0("Número total de genes significativos (FDR < ", input$FDR_MA_T, "): ",
           sum_values(table_MA_T, input$logfc_MA_T, input$FDR_MA_T, input$studies_MA_T)[2])
  )
  
  output$up_genes_MA_T <- renderText(
    paste("Número de genes UP regulados (logFC > ", input$logfc_MA_T, "): ",
          sum_values(table_MA_T, input$logfc_MA_T, input$FDR_MA_T, input$studies_MA_T)[3])
  )
  
  output$down_genes_MA_T <- renderText(
    paste("Numero de genes DOWN regulados (logFC < -", input$logfc_MA_T, "): ",
          sum_values(table_MA_T, input$logfc_MA_T, input$FDR_MA_T, input$studies_MA_T)[4])
  )
  
  
  # Meta-analysis - Blood -----------------------------------------------------
  output$vulcano_MA_B <- renderPlotly(
    vulcano(table_MA_B, input$logfc_MA_B, input$FDR_MA_B, input$studies_MA_B)
  )
  
  output$forest_MA_B <- renderPlotly(
    forest_plot(input$Gene_MA_B, data=table_MA_B, meta_analysis=results_MA_B)
  )
  
  output$table_MA_B <- renderDataTable(
    my_table(table_MA_B, input$logfc_MA_B, input$FDR_MA_B, input$studies_MA_B)
  )
  
  output$download_MA_B <- downloadHandler(
    filename = function(){"metaanalysis.csv"},
    content = function(fname){
      write.csv(my_table(table_MA_B, input$logfc_MA_B, input$FDR_MA_B, input$studies_MA_B), fname)
    }
  )
  
  output$total_genes_MA_B <- renderText(
    paste("Numero total de genes analizados:",
          sum_values(table_MA_B, input$logfc_MA_B, input$FDR_MA_B, input$studies_MA_B)[1])
  )
  
  output$pval_genes_MA_B <- renderText(
    paste0("Número total de genes significativos (FDR < ", input$FDR_MA_B, "): ",
           sum_values(table_MA_B, input$logfc_MA_B, input$FDR_MA_B, input$studies_MA_B)[2])
  )
  
  output$up_genes_MA_B <- renderText(
    paste("Número de genes UP regulados (logFC > ", input$logfc_MA_B, "): ",
          sum_values(table_MA_B, input$logfc_MA_B, input$FDR_MA_B, input$studies_MA_B)[3])
  )
  
  output$down_genes_MA_B <- renderText(
    paste("Numero de genes DOWN regulados (logFC < -", input$logfc_MA_B, "): ",
          sum_values(table_MA_B, input$logfc_MA_B, input$FDR_MA_B, input$studies_MA_B)[4])
  )
  
  # Meta-analysis - Compare ----------------------------------------------------
  output$compare_UP <- renderDataTable(
    compare_table(table_MA_T, table_MA_B, input, first_sign="UP",
                  second_sign="UP")
  )
  
  output$compare_DOWN <- renderDataTable(
    compare_table(table_MA_T, table_MA_B, input, first_sign="DOWN",
                  second_sign="DOWN")
  )
  
  output$compare_UD <- renderDataTable(
    compare_table(table_MA_T, table_MA_B, input, first_sign="UP",
                  second_sign="DOWN")
  )
  
  output$compare_DU <- renderDataTable(
    compare_table(table_MA_T, table_MA_B, input, first_sign="DOWN",
                  second_sign="UP")
  )
}