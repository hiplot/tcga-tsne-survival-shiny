#library(shiny)
library(plotly)
library(shinydashboard)
library(dashboardthemes)

# theme_blue_gradient1 ----------------------------------------------------

theme_blue_gradient1 <- shinyDashboardThemeDIY(
  
  ### general
  appFontFamily = "Arial"
  ,appFontColor = "rgb(0,0,0)"
  ,primaryFontColor = "rgb(0,0,0)"
  ,infoFontColor = "rgb(0,0,0)"
  ,successFontColor = "rgb(0,0,0)"
  ,warningFontColor = "rgb(0,0,0)"
  ,dangerFontColor = "rgb(0,0,0)"
  ,bodyBackColor = "rgb(229,240,253)"
  
  ### header
  ,logoBackColor = "rgb(23,103,124)"
  
  ,headerButtonBackColor = "rgb(238,238,238)"
  ,headerButtonIconColor = "rgb(75,75,75)"
  ,headerButtonBackColorHover = "rgb(210,210,210)"
  ,headerButtonIconColorHover = "rgb(0,0,0)"
  
  ,headerBackColor = "rgb(238,238,238)"
  ,headerBoxShadowColor = "#aaaaaa"
  ,headerBoxShadowSize = "2px 2px 2px"
  
  ### sidebar
  ,sidebarBackColor = cssGradientThreeColors(
    direction = "down"
    ,colorStart = "rgb(20,97,117)"
    ,colorMiddle = "rgb(56,161,187)"
    ,colorEnd = "rgb(3,22,56)"
    ,colorStartPos = 0
    ,colorMiddlePos = 50
    ,colorEndPos = 100
  )
  ,sidebarPadding = 0
  
  ,sidebarMenuBackColor = "transparent"
  ,sidebarMenuPadding = 0
  ,sidebarMenuBorderRadius = 0
  
  ,sidebarShadowRadius = "3px 5px 5px"
  ,sidebarShadowColor = "#aaaaaa"
  
  ,sidebarUserTextColor = "rgb(255,255,255)"
  
  ,sidebarSearchBackColor = "rgb(55,72,80)"
  ,sidebarSearchIconColor = "rgb(153,153,153)"
  ,sidebarSearchBorderColor = "rgb(55,72,80)"
  
  ,sidebarTabTextColor = "rgb(255,255,255)"
  ,sidebarTabTextSize = 13
  ,sidebarTabBorderStyle = "none none solid none"
  ,sidebarTabBorderColor = "rgb(35,106,135)"
  ,sidebarTabBorderWidth = 1
  
  ,sidebarTabBackColorSelected = cssGradientThreeColors(
    direction = "right"
    ,colorStart = "rgba(44,222,235,1)"
    ,colorMiddle = "rgba(44,222,235,1)"
    ,colorEnd = "rgba(0,255,213,1)"
    ,colorStartPos = 0
    ,colorMiddlePos = 30
    ,colorEndPos = 100
  )
  ,sidebarTabTextColorSelected = "rgb(0,0,0)"
  ,sidebarTabRadiusSelected = "0px 20px 20px 0px"
  
  ,sidebarTabBackColorHover = cssGradientThreeColors(
    direction = "right"
    ,colorStart = "rgba(44,222,235,1)"
    ,colorMiddle = "rgba(44,222,235,1)"
    ,colorEnd = "rgba(0,255,213,1)"
    ,colorStartPos = 0
    ,colorMiddlePos = 30
    ,colorEndPos = 100
  )
  ,sidebarTabTextColorHover = "rgb(50,50,50)"
  ,sidebarTabBorderStyleHover = "none none solid none"
  ,sidebarTabBorderColorHover = "rgb(75,126,151)"
  ,sidebarTabBorderWidthHover = 1
  ,sidebarTabRadiusHover = "0px 20px 20px 0px"
  
  ### boxes
  ,boxBackColor = "rgb(255,255,255)"
  ,boxBorderRadius = 5
  ,boxShadowSize = "0px 1px 1px"
  ,boxShadowColor = "rgba(0,0,0,.1)"
  ,boxTitleSize = 16
  ,boxDefaultColor = "rgb(245,133,122)"
  ,boxPrimaryColor = "rgb(115,209,115)"
  ,boxInfoColor = "rgb(242,168,94)"
  ,boxSuccessColor = "rgb(182,129,219)"
  ,boxWarningColor = "rgb(242,127,121)"
  ,boxDangerColor = "rgb(160, 175, 186)"
  
  ,tabBoxTabColor = "rgb(255,255,255)"
  ,tabBoxTabTextSize = 14
  ,tabBoxTabTextColor = "rgb(0,0,0)"
  ,tabBoxTabTextColorSelected = "rgb(0,0,0)"
  ,tabBoxBackColor = "rgb(248,248,248)"
  ,tabBoxHighlightColor = "rgba(44,222,235,1)"
  ,tabBoxBorderRadius = 5
  
  ### inputs
  ,buttonBackColor = "rgb(245,245,245)"
  ,buttonTextColor = "rgb(0,0,0)"
  ,buttonBorderColor = "rgb(200,200,200)"
  ,buttonBorderRadius = 5
  
  ,buttonBackColorHover = "rgb(235,235,235)"
  ,buttonTextColorHover = "rgb(100,100,100)"
  ,buttonBorderColorHover = "rgb(200,200,200)"
  
  ,textboxBackColor = "rgb(255,255,255)"
  ,textboxBorderColor = "rgb(200,200,200)"
  ,textboxBorderRadius = 5
  ,textboxBackColorSelect = "rgb(245,245,245)"
  ,textboxBorderColorSelect = "rgb(200,200,200)"
  
  ### tables
  ,tableBackColor = "rgb(255,255,255"
  ,tableBorderColor = "rgb(240,240,240)"
  ,tableBorderTopSize = 1
  ,tableBorderRowSize = 1
  
)


# ui ----------------------------------------------------------------------
ui <- dashboardPage(
  
  dashboardHeader(title="TCGA Survival Analysis"),
  

  # dashboardSidebar --------------------------------------------------------
  dashboardSidebar(
    sidebarMenu(
      menuItem("About", tabName = "about"),
      menuItem("Pre-generated t-SNE Analysis", tabName = "scatterplots"),
      menuItem("Custom t-SNE / UMAP Analysis" , tabName = "custom")
    )),
  
  # dashboardBody --------------------------------------------------------
  dashboardBody(
    theme_blue_gradient1,
    tags$style(HTML(".main-sidebar { font-size:15px;}")), #change the font size to 20
    tabItems(
      # Tab 1 - About Page  --------------------------------------------------------
      tabItem(tabName = "about",
              column(12, htmlOutput("cancer_name_about")),
              #htmlOutput("inc"),
              includeHTML("www/test.html")
      ),
      # Tab 2 - Original t-SNE Analysis --------------------------------------------------------
      tabItem(tabName = "scatterplots",
              ####### Row 1 - Cancer and Pathway Selector ###########
              fluidRow(
                box(width=12,
                  column(4,uiOutput("cancerSelector")),
                  column(4,conditionalPanel(
                    condition = "typeof input.cancer != 'undefined' && input.cancer != ''",
                    uiOutput("pathway1Selector"))),
                  column(4,conditionalPanel(
                    condition = "typeof input.cancer != 'undefined' && input.cancer != '' && typeof input.pathway1 != 'undefined' && input.pathway1 != ''",
                    uiOutput("pathway2Selector"))),
                  title="Step 1: Select Cancer & Pathway(s) to Analyze",status="danger",solidHeader=TRUE,collapsible=TRUE,background = "light-blue")
              ),
              ####### Row 2 and 3 - Cancer and Pathway Names ###########
              conditionalPanel(condition = "typeof input.cancer != 'undefined'",
                fluidRow(align="center",
                  column(2),
                  column(10,htmlOutput("cancer_name"))),
                # Row 2
                fluidRow(align="center",
                   column(2),
                   column(5, htmlOutput("pathway1_title_3d")),
                   column(5, htmlOutput("pathway2_title_3d"))     
                )
              ),
              ####### Row 4 - 3D Plots and Individial Survival Plots ###########
              conditionalPanel(condition = "typeof input.cancer != 'undefined' && typeof input.pathway1 != 'undefined'",
                fluidRow(align="center",
                      #Step 2 Options (Column 1)   
                      column(2,
                        box(width=12,
                         uiOutput("filter_option"),
                         uiOutput("filter_pheno_selector"),
                          tags$style(type="text/css","
                          .irs-bar {background-color: gray}
                          .irs-bar-edge {background-color: black}
                          .irs-max {color: white}
                          .irs-min {color: white}
                          .irs-single {color:white; background:green;}"),
                         uiOutput("filter_pheno_value1"),
                         uiOutput("filter_comptype_selector1"),
                         uiOutput("filter_pheno_value2"),
                         uiOutput("filter_comptype_selector2")
                        ,title="Step 2: Filter Phenotype(s) (optional)",status="danger",solidHeader=TRUE,collapsible=TRUE,background = "light-blue")
                      ),
                      #Pathway 1 (Column 2)
                      column(5,
                          box(width=NULL,plotlyOutput('scatterplot_out'),tableOutput("counttable_scatter_p1"),title="t-SNE Clustering", status="primary",solidHeader = TRUE,collapsible = TRUE),
                          box(width=NULL,plotOutput('survivalplot_P1_out'),title=" Survival Analysis",status="success",solidHeader = TRUE,collapsible = TRUE)
                      ),
                      #Pathway 2 (Column 3)
                      column(5,
                        conditionalPanel(condition= "typeof input.pathway2 != 'undefined' && input.pathway2 !=''",
                         box(width=NULL,plotlyOutput('scatterplot_out2'),tableOutput("counttable_scatter_p2"),title="t-SNE Clustering", status="primary",solidHeader = TRUE,collapsible = TRUE),
                         box(width=NULL,plotOutput('survivalplot_P2_out'),title=" Survival Analysis",status="success",solidHeader = TRUE,collapsible = TRUE)
              )))),
              ####### Row 5 - Sequential Survival Analysis Using Pathway 1 and 2 ###########
              conditionalPanel(condition = "typeof input.cancer != 'undefined' && typeof input.pathway1 != 'undefined' && typeof input.pathway2 != 'undefined' && input.pathway2 !=''",
                fluidRow(align="center",
                  #Col 1
                  box(width=2,
                     column(width=12,
                      uiOutput("clusterSelector"),
                      uiOutput("clusterSelector2")
                     ),title="Step 3: Sequential Analysis",status="danger",solidHeader=TRUE,collapsible=TRUE,background = "light-blue"
                  ),
                  #Col 2
                  box(width=10,
                    column(12,align="center",
                     htmlOutput("survivalplotheader"),
                     htmlOutput("survivalplotmid"),
                     htmlOutput("survivalplotheader2"),
                     plotOutput('survivalplot_out'),
                     htmlOutput("text2"),
                     #tableOutput("counttable")
                     DT::dataTableOutput("counttable")
                    ),title="Sequential Pathways t-SNE Clusters - Survival Analysis",status="info",solidHeader=TRUE,collapsible=TRUE
              ))),
              ####### Row 6 - Heatmap and Sequential Survival Analysis Using Dendogram ###########
              conditionalPanel(condition = "typeof input.cancer != 'undefined' && typeof input.pathway1 != 'undefined' && input.pathway1 !=''",
                 #1
                 fluidRow(align="center",
                    column(2),
                    box(width=10,
                      column(12,align="center",
                        uiOutput("create_heatmap_button"),
                        plotOutput('heatmap')),
                      title="Heirarchical Clustering of all TCGA Expression Data",status="warning",solidHeader = TRUE,collapsible = TRUE)
                 ),
                 #2
                 fluidRow(align="center",
                    #Col 1
                    box(width=2,
                     column(width = 12,align="center",
                      uiOutput("clusterSelector_heatmap"),
                      uiOutput("groupSelector2")),
                      title="Step 4: Sequential Dendrogram Analysis with Pathway 1",status="danger",solidHeader=TRUE,collapsible=TRUE,background = "light-blue"
                    ),
                    #Col 2
                    box(width=10,
                     column(width = 12,align="center",
                       htmlOutput("heatmaptitle"),
                       htmlOutput("dendrosurvivalplotheader"),
                       htmlOutput("dendrosurvivalplotmid"),
                       htmlOutput("dendrosurvivalplotheader2"),
                       plotOutput('survivalplot_dendro_out'),
                       htmlOutput("text3"),
                       DT::dataTableOutput("counttable_dendro")),
                       title="Sequential Heatmap Dendrogram - Survival Analysis",status="warning",solidHeader = TRUE,collapsible = TRUE)
              )) 
       ),#/tabItem2
      # Tab 3 - Custom t-SNE / UMAP Analysis ------------------------------------
      tabItem(tabName = "custom",
              ####### Row 1 - Cancer and Pathway Selector ###########
              fluidRow(
                box(width=12,
                  column(4,uiOutput("custom_cancerSelector")),
                  column(4,
                    conditionalPanel(condition = "typeof input.custom_cancer != 'undefined' && input.custom_cancer != ''",
                      uiOutput("custom_pathway1_selector"))
                  ),
                  column(4,
                    conditionalPanel(condition = "typeof input.custom_cancer != 'undefined' && input.custom_cancer != '' && typeof input.custom_pathway1 != 'undefined' && input.custom_pathway1 != ''",
                      uiOutput("custom_pathway2_selector"))
                  ),
                  column(4),
                  column(4,conditionalPanel(condition="typeof input.custom_pathway1 != 'undefined' && input.custom_pathway1=='Custom Gene Set'",
                    uiOutput("custom_pathway1_multigene_selector"))
                  ),
                  column(4,conditionalPanel(condition="typeof input.custom_pathway2 != 'undefined' && input.custom_pathway2=='Custom Gene Set'",
                    uiOutput("custom_pathway2_multigene_selector"))
                  ),
                  title="Step 1a: Select Cancer & Pathway(s) to Analyze",status="danger",solidHeader=TRUE,collapsible=TRUE,background = "light-blue")
              ),
              ####### Row 2 - Run Options ###########
              conditionalPanel(condition="typeof input.custom_cancer != 'undefined' && input.custom_cancer != '' && typeof input.custom_pathway1 != 'undefined' && input.custom_pathway1 != ''",
                fluidRow(
                  box(width=12,
                    fluidRow(
                    column(4,uiOutput('exp_data_type')),
                    column(4,uiOutput('run_algo')),
                    ),
                    fluidRow(
                      column(4,
                        conditionalPanel(condition="typeof input.use_UMAP != 'undefined' && input.use_UMAP == 'UMAP'",
                          uiOutput('umap_input_n_epochs')),
                      ),
                      column(4,
                        conditionalPanel(condition="typeof input.use_UMAP != 'undefined' && input.use_UMAP == 'UMAP'",
                          uiOutput('umap_input_min_dist')),
                      ),
                      column(4,
                         conditionalPanel(condition="typeof input.use_UMAP != 'undefined' && input.use_UMAP == 'UMAP'",
                          uiOutput('umap_n_neighbors')),
                    )),
                    fluidRow(
                      column(4,
                       conditionalPanel(condition="typeof input.use_UMAP != 'undefined' && input.use_UMAP == 't-SNE'",
                          uiOutput('tsne_input_max_iter')),
                      ),
                      column(4,
                        conditionalPanel(condition="typeof input.use_UMAP != 'undefined' && input.use_UMAP == 't-SNE'",
                          uiOutput('tsne_input_perplexity')),
                      ),
                      
                    ),
                    title="Step 1b: Select options for dimensionality reduction",status="danger",solidHeader=TRUE,collapsible=TRUE,background = "light-blue")
                ),
              ####### Row 3 and 4 - Cancer and Pathway Names ###########
                fluidRow(align="center",
                  column(2),
                  column(10,htmlOutput("custom_cancer_name"))
                ),
                fluidRow(align="center",
                 column(2),
                 column(5, htmlOutput("custom_pathway1_title_3d")),
                 column(5, htmlOutput("custom_pathway2_title_3d"))     
                ),
              ####### Row 5 - 3D Plots and Individial Survival Plots ###########
                fluidRow(align="center",
                  #Step 2: Options (Column 1)   
                  column(2,
                   box(width=12,
                    htmlOutput("text_filter_custom"),
                    uiOutput("filter_option_custom"),
                    uiOutput("filter_pheno_selector_custom"),
                    tags$style(type="text/css","
                    .irs-bar {background-color: gray}
                    .irs-bar-edge {background-color: black}
                    .irs-max {color: white}
                    .irs-min {color: white}
                    .irs-single {color:white; background:green;}"),
                    uiOutput("filter_custom_pheno_value1"),
                    uiOutput("filter_comptype_custom_selector1"),
                    htmlOutput("text_run_custom_warning"),
                    uiOutput("run_custom"),
                    title="Step 2: Filter Phenotype(s) (Optional) and Run Analysis",status="danger",solidHeader=TRUE,collapsible=TRUE,background = "light-blue")
                  ),
                  column(5,
                   box(width=NULL,plotlyOutput('custom_scatterplot_P1_out'),tableOutput("counttable_custom_scatter_p1"),title="Clustering: Pathway 1", status="primary",solidHeader = TRUE,collapsible = TRUE),
                   box(width=NULL,plotOutput('custom_survivalplot_P1_out'),title=" Survival Analysis: Pathway 1",status="success",solidHeader = TRUE,collapsible = TRUE)
                  ),   
                  column(5,
                    conditionalPanel(condition= "typeof input.custom_pathway2 != 'undefined' && input.custom_pathway2 !=''",
                      box(width=NULL,plotlyOutput('custom_scatterplot_P2_out'),tableOutput("counttable_custom_scatter_p2"),title="Clustering: Pathway 2", status="primary",solidHeader = TRUE,collapsible = TRUE),
                      box(width=NULL,plotOutput('custom_survivalplot_P2_out'),title=" Survival Analysis: Pathway 2",status="success",solidHeader = TRUE,collapsible = TRUE)
               )))),
              ####### Row 6 - Sequential Survival Analysis Using Pathway 1 and 2 ###########
              conditionalPanel(condition = "typeof input.custom_cancer != 'undefined' && typeof input.custom_pathway1 != 'undefined' && typeof input.custom_pathway2 != 'undefined' && input.custom_pathway2 !=''",
                fluidRow(align="center",
                  # Col 1
                   box(width=2,
                    column(12,
                      uiOutput("custom_group_clusterSelector"),
                      uiOutput("custom_group_clusterSelector2")
                    ),title="Step 3: Sequential Analysis",status="danger",solidHeader=TRUE,collapsible=TRUE,background = "light-blue"
                  ),
                  # Col 2
                  box(width=10,
                    column(10,align="center",       
                      plotOutput('custom_survivalplot_out'), 
                      DT::dataTableOutput("custom_counttable")
                    ),title="Sequential Pathway Surivial Analysis", status="info",solidHeader = TRUE,collapsible = TRUE)
              )),
              ####### Row 7 - Heatmap and Sequential Survival Analysis Using Dendogram ###########
              conditionalPanel(condition= "typeof input.custom_cancer != 'undefined' && typeof input.custom_pathway1 != 'undefined' && input.custom_pathway1 !=''",
                #1
                fluidRow(align="center",
                  column(2),
                  box(width=10,
                    column(12,align="center",
                      uiOutput("custom_create_heatmap_button"),
                      plotOutput('custom_heatmap')),
                      title="Heirarchical Clustering of all TCGA Expression Data",status="warning",solidHeader = TRUE,collapsible = TRUE)
                ),
                #2
                fluidRow(align="center",
                  #Col 1
                  box(width=2,
                   column(width = 12,align="center",
                      uiOutput("custom_clusterSelector_heatmap"),
                      uiOutput("custom_groupSelector2")),
                    title="Step 4: Sequential Dendrogram Analysis with Pathway 1",status="danger",solidHeader=TRUE,collapsible=TRUE,background = "light-blue"
                  ),
                  #Col 2
                  box(width=10,
                    column(width = 12,align="center",
                    #htmlOutput("heatmaptitle"),
                    #htmlOutput("dendrosurvivalplotheader"),
                    #htmlOutput("dendrosurvivalplotmid"),
                    #htmlOutput("dendrosurvivalplotheader2"),
                    plotOutput('custom_survivalplot_dendro_out'),
                    #htmlOutput("text3"),
                    DT::dataTableOutput("custom_counttable_dendro")),
                    title="Sequential Heatmap Dendrogram - Survival Analysis",status="warning",solidHeader = TRUE,collapsible = TRUE)
              ))
      )#/tabItem 3
  
# closing brackets --------------------------------------------------------
    )#tabItems
  ) #dashboardBody
)#dashboardPage
