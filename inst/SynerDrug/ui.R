require(shiny)
require(shinydashboard)
require(ggplot2)


header <- shinydashboard::dashboardHeader(title = "SynerDrug")

sidebar <- shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
        menuItem("Input", tabName = "indata", icon = icon("file")),
        menuItem("Data", tabName = "data", icon = icon("table")),
        menuItem("Individual responses", tabName = "indResp", icon = icon("line-chart")),
        menuItem("HSA", tabName = "hsa", icon = icon("table")),
        menuItem("Bliss", tabName = "bliss", icon = icon("table")),
        menuItem("Chou-Talalay", tabName = "chou", icon = icon("table")),
        menuItem("Loewe", tabName = "loewe", icon = icon("table")),
        menuItem("Combination/effect", tabName = "combInd", icon = icon("map")),
        menuItem("Global report", tabName = "report", icon = icon("file")) #,
        #menuItem("Help", tabName = "help", icon = icon("question"))
    )
)

body <- shinydashboard::dashboardBody(
    withMathJax(),
    tabItems(
        ##------------------------------------------------------------
        ## data input
        tabItem(tabName = "indata", textOutput("drcOK"), h2("Data input"), br(), #actionButton("dataFormatHelp", "Data Format"),
                fluidRow(
                    column(width = 4, selectInput("inputFormat", "Data input format", choices = c("Matrix with dose" = "matdose", "Long format" = "long"))),
                    column(width = 4, selectInput("content", "Data content", choices = c("Death rate" = "death", "Survival rate" = "surv")))),
                fluidRow(
                    column(width = 4, selectInput("dataRange", "Data range", choices = c("0-1" = "1", "0-100" = "100"))),
                    column(width = 4, selectInput("unit", "Dose unit", choices = c("mM" = "milli", "μM" = "micro", "nM" = "nano")))
                    ),
                fluidRow(
                    column(width = 4, textInput("dA", "Drug A name", "")),
                    column(width = 4, textInput("dB", "Drug B name", ""))
                    ),
                fluidRow(column(width = 4, fileInput("data", label = "File"))),
                uiOutput("dataInfos")
                ),

        ##------------------------------------------------------------
        ## data vis
        tabItem(tabName = "data", h2("Data vis"),
                tabsetPanel(
                    tabPanel("Table", br(), div(DT::dataTableOutput("rawData")), style='width:800px'),
                    tabPanel("Heatmap", #checkboxInput("mergedHeat", "Merge replicates", TRUE),
                             plotOutput("heatmap", width="auto", height="auto")),
                    tabPanel("Parallel plot", uiOutput("parplotSelect"), plotOutput("parplot")),
                    tabPanel("Surface",
                             fluidRow(column(width=4, sliderInput("theta", "azimuth", min=-180, max=180, value=-100, step=5)),
                             column(width=4, sliderInput("phi", "colatitude", min=0, max=50, value=10, step=1))),
                    fluidRow(plotOutput("surfPlot")))
                    ##   tabPanel("Predict effect", fluidRow(column(width=4, uiOutput("userDoseAUI"), uiOutput("userDoseBUI")), column(width=4, br(), br(), htmlOutput("predEffect"))), plotOutput("predEffectPlot"))
                    )
                ),

        ##------------------------------------------------------------
        ## ind resp
        tabItem(tabName = "indResp", h2("Individual responses"),
                tabsetPanel(
                    tabPanel("Curves", plotOutput("plotA"), plotOutput("plotB")),
                    tabPanel("Hill model parameters", br(),# actionButton("hilleq", "Hill Model"),
                             br(), div(DT::dataTableOutput("paramHills")), style='width:100px')#,
                    #tabPanel("Dose calculator", br(),
                     #        fluidRow( column(width=4, numericInput("doseCalc", "Expected effect:", value = 0.5)), column(width=4,tableOutput("EDout"))),
                      #       plotOutput("EDplot"))
                    )
                ),
        ##------------------------------------------------------------
        ## Loewe
        tabItem(tabName = "loewe", h2("Loewe model"), tabsetPanel(
                tabPanel("Data under Loewe model", plotOutput("Loewe")),
                tabPanel("Loewe excess", plotOutput("LoeweExcess")),
                tabPanel("Isobologram",  fluidRow(column(width = 4, numericInput("effectIso", "Effect", value = 0.5, min = 0, max = 1, step = 0.1)), column(width = 4, selectInput("isobolType", "Isobologram type", choices = c("Linear" = "linear", "Steel and Peckham" = "SP")))), plotOutput("isobole"))
            )),

        ##------------------------------------------------------------
        tabItem(tabName = "hsa", h2("Highest single agent"), plotOutput("HSA")),

        ##------------------------------------------------------------
        tabItem(tabName = "bliss", h2("Bliss model"), plotOutput("Bliss")),

        ##------------------------------------------------------------
        tabItem(tabName = "chou", h2("Chou-Talalay"), plotOutput("Chou")),

        ##------------------------------------------------------------
        tabItem(tabName = "combInd", selectInput("combMeth", "Combination Index:", choices = c("HSA" = "HSA", "Bliss" = "Bliss", "Loewe excess" = "LoeweExcess", "Chou-Talalay" = "Chou")), plotOutput("combEffect", brush = "EDbrush"), br(), downloadButton('tableDL', 'Download table'), br(), DT::dataTableOutput("dataCombEffect")),

        ##------------------------------------------------------------
        tabItem(tabName = "report", h2("Results report"), br(),
                checkboxGroupInput(inputId = "reportChoices", label = "Choices", choices = c("Data" = "Data", "HSA" , "Bliss" = "Bliss","Individual responses" = "indResp", "Loewe predicted" = "LoewePred", "Loewe excess" = "LoeweExcess", "Chou-Talalay" = "Chou", "Isobologram" = "isobol"), selected = c("Data", "HSA", "Bliss", "indResp", "LoewePred", "LoeweExcess", "Chou", "isobol"), inline = TRUE),
                uiOutput("isobolUI"),
                fluidRow(column(width = 4, uiOutput("reportTitleUI"))),#, column(width = 4, "BB")),
                br(),
                downloadButton("downloadReport", "Download")),

        ##------------------------------------------------------------
        tabItem(tabName = "help", uiOutput("help"))

        )
    )



dashboardPage(
    header, sidebar, body
)

