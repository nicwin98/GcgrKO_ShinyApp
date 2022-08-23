library(EnhancedVolcano)
library(tidyverse)
library(DT)
library(shiny)
library(shinyjs)
library(shinythemes) 
library(writexl)

KO_WT_shiny <- read.csv("KO_WT_shiny.csv")
t_box_shiny <- read.csv("t_box_shiny.csv") 
GOBP_KO <- read.csv("GOBP_Genes_KOvsWT.csv")
GOBP_KO$Treatment <- rep("GcgrKO", length(GOBP_KO$Description))
GOBP_KO$GeneFraction <- sapply(GOBP_KO$GeneRatio, function(x) eval(parse(text=x)))
t_box_shiny$Condition <- factor(t_box_shiny$Condition, levels = c("GcgrKO","WT"))

ui <- fluidPage(theme = shinytheme("paper"),
            navbarPage(title = "Differential expression analysis - GcgrKO",       
                       tabPanel(title = "Preface",
                                fluidRow(column(12, wellPanel(tags$h4("Transcriptomic analysis of livers from glucagon receptor knock-out (GcgrKO) mice"),
                                                              tags$br(),
                                                              tags$h6("This app contains information on RNA sequencing data from livers of GcgrKO mice compared to wildtype littermates. 
                                                                      For detailed information on the bioinformatic analysis used to generate this data please refer to the R codes 
                                                                      available at", tags$a("https://github.com/nicwin98/GCGA_GCGR-Ab_GcgrKO"), ".
                                                                      Raw data are available through ArrayExpress via the accession code", "and the publication listed below:"),
                                                              tags$br(),
                                                              tags$i(tags$h5("Glucagon acutely regulates hepatic amino acid catabolism and the effect may be disturbed by steatosis")),
                                                              tags$h6("Marie Winther-Sørensen, Katrine D. Galsgaard, Alberto Santos, Samuel A.J. Trammell, Karolina Sulek, 
                                                                      Rune E. Kuhre, Jens Pedersen, Daniel B. Andersen, Anna S. Hassing, Morten Dall, Jonas T. Treebak, 
                                                                      Matthew P. Gillum, Signe S. Torekov, Johanne A. Windeløv, Jenna E. Hunt, Sasha A.S. Kjeldsen, 
                                                                      Sara L. Jepsen, Catherine G. Vasilopoulou, Filip K. Knop, Cathrine Ørskov, Mikkel P. Werge, 
                                                                      Hanne Cathrine Bisgaard, Peter Lykke Eriksen, Hendrik Vilstrup, Lise Lotte Gluud, Jens J. Holst, 
                                                                      and Nicolai J. Wewer Albrechtsen"),
                                                              tags$h6("Molecular Metabolism (2020), DOI:", tags$a("https://doi.org/10.1016/j.molmet.2020.101080")),
                                                              tags$br(),
                                                              tags$h6("Raw data are available through ArrayExpress via the accession code", tags$b("E-MTAB-12060"), ".")))),
                                fluidRow(column(12, wellPanel(
                                  tags$h5("The app contains three tabs where you can find the following information:"),
                                  tags$ul(
                                    tags$li(tags$strong("GcgrKO vs WT:")),
                                    tags$ul(
                                      tags$li("An adjustable volcano plot of all analyzed genes."),
                                      tags$li("A table of all genes displaying the results of the differential expression analysis."),
                                      tags$ul(
                                        tags$li("The information in the table is available for download.")))),
                                  tags$ul(
                                    tags$li(tags$strong("Single Gene Expression:")),
                                    tags$ul(
                                      tags$li("Boxplot showing the sample variation for individual genes."),
                                      tags$ul(
                                        tags$li("The boxplot and datapoints are available for download.")))),
                                  tags$ul(
                                    tags$li(tags$strong("Gene Ontology Biological Pathways (GOBPs)")),
                                    tags$ul(
                                      tags$li("All significant GOBPs are displayed in a table."),
                                      tags$li("Create a dotplot by clicking on your GOBPs of interest."),
                                      tags$ul(
                                        tags$li("The table information and the dotplot are available for download.")))),
                                  tags$br(),
                                  tags$h6("Any inquires related to the application and publication should be directed 
                                                              to the corresponding author, Nicolai J. Wewer Albrechtsen: nicolai.albrechtsen@sund.ku.dk."),
                                  tags$h6("Two other apps have been created in relation to the abovementioned publication:", tags$a("https://weweralbrechtsenlab.shinyapps.io/GCGA"), "&", 
                                          tags$a("https://weweralbrechtsenlab.shinyapps.io/GCGR_Ab"), "."))))),
                       
                       tabPanel(title = "GcgrKO vs WT",
                               fluidRow(column(4,
                                wellPanel(
                                sliderInput(
                                 inputId = "num1",
                                 label = "Choose a log2FoldChange cutoff",
                                 value = 1,
                                 min = 0.0,
                                 max = 5.0,
                                 step = 0.5),
                                sliderInput(
                                  inputId = "xlim1",
                                  label = "Adjust x-axis",
                                  value = c(-24,26),
                                  min = -24,
                                  max = 26,
                                  step = 1),
                                selectInput(
                                  inputId = "radio1",
                                  label = "Choose FDR adjusted p-value:",
                                  selected = 0.05,
                                  choices = list("0.05","0.01","0.001","0.0001","0.00001")))),
                                column(8,plotOutput("volcano"))),
                               fluidRow(column(12, wellPanel(tags$h4("Use the table below to search for genes of interest"),
                                                             tags$h5("You can download the information in the table for selected genes by clicking on your genes of interest and clicking 
                                                                 the download button below the table"),
                                                             tags$h6("The table displays information on:"),
                                                             tags$ul(
                                                               tags$li("ENSEMBL ID: Unique identifier (used to select genes in the Single Gene Expression tab)"),
                                                               tags$li("Gene Symbol"),
                                                               tags$li("Gene Name"),
                                                               tags$li("Base Mean: Mean of normalized expression (not comparable across genes)"),
                                                               tags$li("Log 2 fold changes (green: > 0 | red: < 0)"),
                                                               tags$li("Fold changes (green: > 1 | red: < 1)"),
                                                               tags$li("lfcSE: Standard error of the log 2 fold change"),
                                                               tags$li("Stat: Wald statistic used to calculate the p-values"),
                                                               tags$li("P-value: Not corrected for multiple testing"),
                                                               tags$li("Padj: FDR adjusted p-value (green: padj < 0.05 | black: padj > 0.05)")),
                                                             tags$h6("Sometimes adjusted p-values or p-values and adjusted p-values are missing. This is likely due to:"),
                                                             tags$ul(
                                                               tags$li("Missing padj: The gene is to lowly expressed and thus not tested (see baseMean). 
                                                                   DESeq2's independent filtering function removed it."),
                                                               tags$li("Missing pvalue and padj: The gene has one or more extreme outliers. 
                                                                   As a result the gene is removed from the stastical analysis
                                                                   (see if this is the case by using the Single Gene Expression tab).")),
                                                             tags$h6("Green log 2 fold changes / fold changes means the gene is up-regulated in the GcgrKO mice compared
                                                                     to the wildtype littermates.")))),
                               fluidRow(column(12, wellPanel(DTOutput("listDE")))),
                               fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsDE",
                                                                          label = "Clear selected rows"),
                                                             downloadButton("printDEselect", 'Download info on all selected genes'))))),
                  
                  tabPanel(title = "Single Gene Expression - BoxPlot",
                           fluidRow(column(12, wellPanel(tags$h5("The boxplot below displays expression counts normalized with DESeq2 (plotCount function) across treatment groups"),
                                                         tags$h6("The counts are not comparable across genes, as they are not normalized to gene length. The boxplot shows the median, 25", 
                                                                 tags$sup("th"),", and 75", tags$sup("th"),"percentiles. Points are displayed as outliers if they are above or below 
                                                                 1.5 times the interquartile range.")))),
                           fluidRow(column(12,
                                           wellPanel(
                                             textInput(
                                               inputId = "num5",
                                               label = "Find the ENSEMBL ID of your gene (see the table under GcgrKO vs WT), enter it below, and click the update button.
                                                       Remember to change the title accordingly.",
                                               value = "ENSMUSG00000025991.9"),
                                             textInput(
                                               inputId = "plotboxtitle",
                                               label = "Write your desired title",
                                               value = "CPS1 expression - ENSMUSG00000025991.9"),
                                             actionButton(
                                               inputId = "updateBox",
                                               label = "Update"),
                                             downloadButton("printboxplot", 'Download plot as a PDF file'),
                                             downloadButton("printdatatable", 'Download expression values from plot')))),
                           fluidRow(column(12, plotOutput("expressionboxplot")))),
                  tabPanel(title = "Gene Ontology Biological Process",
                           fluidRow(column(12, wellPanel(h5("The table below allows you to click and select the Gene Ontology Biological Pathways (GOBPs) you are interested in. 
                                                            The selected GOBPs will appear in the dotplot below the table."),
                                                         h6("Both the displayed dotplot and the selected table information are available for download. 
                                                            All GOBPs in the table are significantly enriched in the treatment group compared to the control group.")))),
                           fluidRow(column(4,
                                           downloadButton("printGOBPselect", 'Download info on all selected GOBPs'),
                                           downloadButton("printGOBPdotplot", "Download the current dotplot")),
                                    column(4,
                                           numericInput(
                                             inputId = "radio2",
                                             label = "Adjust the width of the downloadable dotplot",
                                             value = 8,
                                             min = 1,
                                             max = 15,
                                             width = "200%")),
                                    column(4,
                                           numericInput(
                                             inputId = "radio3",
                                             label = "Adjust the height of the downloadable dotplot",
                                             value = 4,
                                             min = 1,
                                             max = 15,
                                             width = "200%"))),
                           fluidRow(column(12, wellPanel(DTOutput("GOBPtable")))),
                           fluidRow(column(12, wellPanel(actionButton(inputId = "clearRowsGOBP",
                                                                      label = "Clear selected rows")))),
                           fluidRow(useShinyjs(),style = "padding-top:20px", column(12, plotOutput("GOBPlot"))))
                  )
            )

server <- function(input, output, server) {
  
  data_box <- eventReactive(input$updateBox, {filter(t_box_shiny, ENSEMBL == input$num5, .preserve = TRUE)})
  
  reactive_KO_WT <- reactive({KO_WT_shiny})
  
  output$volcano <- renderPlot({
    EnhancedVolcano(reactive_KO_WT(),
                    lab = KO_WT_shiny$gene_symbol,
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = 'GcgrKO vs WT',
                    pCutoff = as.numeric(input$radio1),
                    FCcutoff = input$num1,
                    pointSize = 1.0,
                    labSize = 2.0,
                    cutoffLineWidth = 0.5,
                    col=c('darkgray', 'darkgray', 'blue', 'red'),
                    colAlpha = 0.75,
                    legendLabels=c('Not sig.','','Sig. & low log2FC','Sig. & high log2FC'),
                    drawConnectors = TRUE,
                    widthConnectors = 0.5,
                    xlim = as.numeric(input$xlim1),
                    ylim = c(0,118))
  })
  
  output$listDE <- renderDT({datatable(KO_WT_shiny, 
                                       options = list(
                                         scrollX = TRUE,
                                         pageLength = 10,
                                         lengthMenu = c(5,10,25,50,200),
                                         filter = "bottom"
                                       )) %>%
      formatSignif(., columns = c(9,10), digits = 4) %>%
      formatRound(., columns = c(4:8), digits = 4, interval = 0, mark = "") %>%
      formatStyle('padj', color = styleInterval(0.05,c('green','black'))) %>%
      formatStyle('log2FoldChange', color = styleInterval(0, c('red', 'green'))) %>%
      formatStyle('foldChange', color = styleInterval(1, c('red', 'green')))
  })
  
  DTreset <- DT::dataTableProxy("listDE")
  shiny::observeEvent(input$clearRowsDE, {
    DT::selectRows(DTreset, NULL)
  })
  
  output$printDEselect = downloadHandler('Gcgr_KO_Genes_Selected.xlsx', content = function(file) 
  {srows_data <- input$listDE_rows_selected
  writexl::write_xlsx(KO_WT_shiny[srows_data, , drop = FALSE], path = file)
  })
  
  output$expressionboxplot <- renderPlot({
    ebplot <- ggplot(data_box(),aes(x=Condition, y=count, fill = Condition)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(col = Condition), alpha = 0.6, color = "black", width = 0.15) +
      ylab("Count") +
      xlab("Condition") +
      ggtitle(input$plotboxtitle) +
      scale_fill_manual(values = c("#FF3333","#CCCCCC")) +
      theme_bw()
    print(ebplot)
  })
  
  output$printboxplot <- downloadHandler(
    filename = function() {
      paste('ExpressionBoxPlot_GCGR_Analog', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, 
             ggplot(data_box(),aes(x=Condition, y=count, fill = Condition)) +
               geom_boxplot(outlier.shape = NA) +
               geom_jitter(aes(col = Condition), alpha = 0.6, color = "black", width = 0.15) +
               ylab("Count") +
               xlab("Condition") +
               ggtitle(input$plotboxtitle) +
               scale_fill_manual(values = c("#FF3333","#CCCCCC")) +
               theme_bw() , 
             dpi = 300)
    })
  
  output$printdatatable <- downloadHandler(
    filename = function() {paste('NormCounts_GcgrKO_', Sys.Date(), '.xlsx', sep='')},
    content = function(file) 
      {writexl::write_xlsx(dplyr::filter(data_box()), path = file)
    })
  
  output$GOBPtable <- renderDT(GOBP_KO,
                               filter = "bottom",
                               options = list(
                                 autoWidth = TRUE,
                                 scrollX = TRUE,
                                 rowCallback = JS(
                                   "function(row, data) {",
                                   "for (i = 1; i < data.length; i++) {",
                                   "if (data[i]<0.01){",
                                   "$('td:eq('+i+')', row).html(data[i].toExponential(4));",
                                   "}",
                                   "}",
                                   "}"),
                                 columnDefs = list(list(width = '275px', targets = c(3)),
                                                   list(width = '100px', targets = c(1)),
                                                   list(width = '55px', targets = c(6,7)),
                                                   list(visible=FALSE, targets=c(9,10))),
                                 pageLength = 10,
                                 lengthMenu = c(5,10,25,50,200))
  )
  
  DTreset2 <- DT::dataTableProxy("GOBPtable")
  shiny::observeEvent(input$clearRowsGOBP, {
    DT::selectRows(DTreset2, NULL)
  })
  
  output$GOBPlot <- renderPlot({
    validate(
      need(input$GOBPtable_rows_selected, "Please click and select your GOBPs of interest to generate the dotplot.")
    )
    srows_data <- input$GOBPtable_rows_selected
    GOBP_KO <- GOBP_KO[srows_data, , drop = FALSE]
    ggplot(GOBP_KO, aes(x = Genes, y = Description)) + 
      geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
      scale_color_distiller(palette = "Spectral", direction = 1) +
      ylab(NULL) +
      xlab(NULL) +
      facet_wrap(~Treatment, drop = TRUE) +
      theme_bw(base_size = 14) +
      theme(axis.text = element_text(color = "black", size = 12), 
            legend.key.size = unit(0.7, 'cm'), 
            legend.title = element_text(size=11), 
            legend.text = element_text(size=9)) +
      guides(size = guide_legend(order=1))
  })
  
  output$printGOBPselect = downloadHandler('GcgrKO_GOBP_Selected.xlsx', content = function(file) 
  {srows_data <- input$GOBPtable_rows_selected
  writexl::write_xlsx(GOBP_KO[srows_data, , drop = FALSE], path = file)
  })
  
  output$printGOBPdotplot <- downloadHandler(
    filename = function() {
      paste('DotPlot_GcgrKO_GOBP', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      srows_data <- input$GOBPtable_rows_selected
      GOBP_KO <- GOBP_KO[srows_data, , drop = FALSE]
      ggsave(file, 
             ggplot(GOBP_KO, aes(x = Genes, y = Description)) + 
               geom_point(aes(size = GeneFraction, color = p.adjust), stroke = 2) +
               scale_color_distiller(palette = "Spectral", direction = 1) +
               ylab(NULL) +
               xlab(NULL) +
               facet_wrap(~Treatment, drop = TRUE) +
               theme_bw(base_size = 14) +
               theme(axis.text = element_text(color = "black", size = 12), 
                     legend.key.size = unit(0.7, 'cm'), 
                     legend.title = element_text(size=11), 
                     legend.text = element_text(size=9)) +
               guides(size = guide_legend(order=1)), 
             dpi = 1000, width = as.numeric(input$radio2), height = as.numeric(input$radio3))
    })
  
}

shinyApp(ui, server)