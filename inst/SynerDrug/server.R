shinyServer(function(input, output){

    ########################################
    ## Input
    ########################################
    datasetInput <- reactive({
        inFile <- input$data
        if (input$inputFormat == "long") {
            header <- TRUE
        } else header <- FALSE
        if (!is.null(inFile)) {
            dec <- synerdrug:::detectDec(inFile$datapath)
            D <- read.table(inFile$datapath, sep = "\t", header = header, dec = dec)
        }
        else D <- NULL
        return(D)
    })

    dataUnit <- reactive({
        return(switch(input$unit, milli = "(mM)", micro = "(μM)", nano = "(nM)"))
    })

    ## Create DrugSyn object
    dataLong <- reactive({
        D <- datasetInput()
        if (!is.null(D)) {
            if (input$inputFormat == "matdose"){
                NAlines <-  which(apply(D, 1, function(x) all(is.na(x))))
                if (length(NAlines) == 0) NAlines <- nrow(D) + 1
                nB <- ncol(D) - 1
                nA <- NAlines[1] - 2
                doseB <- unlist(D[1,-1])
                D <- D[-c(1, NAlines, NAlines+1), ]
                doseA <- D[,1]
                nRep <- nrow(D)/nA
                D <- reshape2::melt(D, id.vars="V1")
                D$concB <- doseB[D$variable]
                D$variable <- NULL
                D$rep <- rep(rep(1:nRep, each = nA), nB)
                D$concA <- D$V1
                D$V1 <- NULL
                D$PosA <- match(D$concA, sort(unique(D$concA)))
                D$PosB <- match(D$concB, sort(unique(D$concB)))
                D$A <- D$PosA
                D$B <- D$PosB
            colnames(D) <- c("value", "concB", "rep", "concA", "PosA", "PosB", drugA(), drugB())
            } else {
                D$PosA <- match(D[, drugA()], sort(unique(D[, drugA()])))
                D$PosB <- match(D[, drugB()], sort(unique(D[, drugB()])))
                D$concA <- D[, drugA()]
                D$concB <- D[, drugB()]
                doseA <- unique(D$concA)
                doseB <- unique(D$concB)
                D$concAOrder <- match(D$concA, sort(doseA))
                D$concBOrder <- match(D$concB, sort(doseB))
                D[, drugA()] <- D$PosA
                D[, drugB()] <- D$PosB
            }
            D <- D[, c("value", "rep", "concA", "concB", "PosA", "PosB", drugA(), drugB())]
            D$valueOrig <- D$value
            if (input$dataRange == 100) D$value <- D$value/100
            colnames(D) <- c("value", "rep", paste0("conc", drugA()), paste0("conc", drugB()), paste0("Pos", drugA()), paste0("Pos", drugA()), drugA(), drugB(), "valueOrig")
        }
        return(D)
    })

    drugA <- reactive(input$dA)
    drugB <- reactive(input$dB)

    content <- reactive(switch(input$content, death = "Death", surv = "Survival"))

    object <- reactive({
        D <- dataLong()
        if (!is.null(D)) {
            doses <- list(a = sort(unique(D[, paste0("conc", drugA())])), b = sort(unique(D[, paste0("conc", drugB())])))
            names(doses) <- c(drugA(), drugB())
            D <- D[, c("value", "rep", paste0("conc", c(drugA(), drugB())))]
            colnames(D)[3:4] <- c(drugA(), drugB())
            object <- synerdrug:::makeDrugSyn(D, doses = doses, content = content(), typeHill = 4)
        } else object <- NULL
        return(object)
    })


    #########################################
    ## retrieve design infos
    #########################################
    ## doses
    doseA <- reactive({
        D <- object()
        if (is.null(D)) return(NULL)
        else return(doses(D)[[1]])
    })
    doseB <- reactive({
        D <- object()
        if (is.null(D)) return(NULL)
        else return(doses(D)[[2]])
    })
    ## doses number
    nA <- reactive({
        D <- object()
        if (is.null(D)) return(NULL)
        else return(length(doseA()))
    })
    nB <- reactive({
        D <- object()
        if (is.null(D)) return(NULL)
        else return(length(doseB()))
    })

    ## dose information on Input page
    output$doseUIA <- renderUI({
        if (is.null(nA())) return(NULL)
        else return(lapply(1:nA(), function(x) textInput(paste0("doseA",x), label = paste("Dose", x), value = doseA()[x])))
    })

    output$doseUIB <- renderUI({
        if (is.null(nB())) return(NULL)
        else return(lapply(1:nB(), function(x) textInput(paste0("doseB",x), label = paste("Dose", x), value = doseB()[x])))
    })

    output$dataInfos <- renderUI({
        if (is.null(object())) return(NULL)
        else {
            res <- list(
                fluidRow(
                    column(width = 4, uiOutput("doseUIA")),
                    column(width = 4, uiOutput("doseUIB")))
            )
            return(res)
        }
    })


    #####################################
    ## Data
    #####################################
    ## data table visualisation
    output$rawData <- DT::renderDataTable({
        d <- round(dataLong(), digits = 3)
        d <- d[, c(drugA(), drugB(), paste0("conc", c(drugA(), drugB())), "rep", "value")]
        return(DT::datatable(d, rownames = FALSE, options = list(dom="lrtip")))
    })

    ## heatmap
    output$heatmap <- renderPlot({
        ## TODO rep séparés
        ## if (input$mergedHeat)
        g <- plot(object(), type = "heatmap", what = "value")
        print(g)
    }, bg="transparent", width=600, height=600)

    ## select input for parallel plot
    output$parplotSelect <- renderUI({
        choices <- list(1, 2)
        names(choices) <- c(drugA(), drugB())
        selectInput("parplotWay", label = "select", choices = choices)
    })

    ## Parallel lines plot
    output$parplot <- renderPlot({
        ref <- c(drugA(), drugB())[as.numeric(input$parplotWay)]
        g <- plot(object(), type = "parallel", ref = ref)
        g <- g + ggplot2::theme(plot.background = ggplot2::element_blank(), legend.background = ggplot2::element_blank())
        print(g)
    }, bg = "transparent")

    ## Surface
    output$surfPlot <- renderPlot({
        plot(object(), type = "surface", what = "value",  phi = input$phi, theta = input$theta, expand = 0.8)
    }, bg = "transparent", height = 600, width = 600)


    #################################
    ## Individual response
    #################################
    output$plotA <- renderPlot({
        g <- plot(object(), type = "ind", drug = drugA())
        g <- g + ggplot2::theme(plot.background = ggplot2::element_blank(), legend.background = ggplot2::element_blank())
        print(g)
    }, bg = "transparent")
    output$plotB <- renderPlot({
        g <- plot(object(), type = "ind", drug = drugB())
        g <- g + ggplot2::theme(plot.background = ggplot2::element_blank(), legend.background = ggplot2::element_blank())
        print(g)
    }, bg = "transparent")

    output$paramHills <- DT::renderDataTable({
        d <- as.data.frame(respInd(object()))
        return(DT::datatable(d, rownames = TRUE, options = list(dom = "", autoWidth = TRUE, columnDefs = list(list(width = '30px', targets = 1:2, sort=FALSE)))))
    })

    #################################
    ## HSA
    #################################
    output$HSA <- renderPlot({
        g <- plot(object(), type = "heatmap", what = "HSA")
        print(g)
    }, bg = "transparent", width = 600, height = 600)

    #################################
    ## Bliss
    #################################
    output$Bliss <- renderPlot({
        g <- plot(object(), type = "heatmap", what = "Bliss")
        print(g)
    }, bg = "transparent", width = 600, height = 600)

    #################################
    ## Chou
    #################################
    output$Chou <- renderPlot({
        g <- plot(object(), type = "heatmap", what = "Chou")
        print(g)
    }, bg = "transparent", width = 600, height = 600)

    ################################
    ## Loewe
    ################################
    output$Loewe <- renderPlot({
        g <- plot(object(), type = "heatmap", what = "Loewe")
        print(g)
    }, bg = "transparent", width = 600, height = 600)
    output$LoeweExcess <- renderPlot({
        g <- plot(object(), type = "heatmap", what = "LoeweExcess")
        print(g)
    }, bg = "transparent", width = 600, height = 600)

    output$isobole <- renderPlot({
        g <- isobologram(object(), effect = input$effectIso, mode = input$isobolType)
        g <- g + ggplot2::theme(plot.background = ggplot2::element_blank(), legend.background = ggplot2::element_blank())
        print(g)
    }, bg = "transparent", width = 600, height = 600)

    ################################
    ## Combination / effect
    ################################
    ## data for plot
    dataEDplot <- reactive({
        d <- synerdrug:::meanData(object())
        d <- d[d[, drugA()] != 0 & d[, drugB()] != 0, ] ## filter mono
        rownames(d) <- as.character(1:nrow(d))
        return(d)
    })
    dataEDplotBrushed <- reactive({
        d <- dataEDplot()
        if (!is.null(input$EDbrush) && nrow(brushedPoints(d, input$EDbrush))) d <- brushedPoints(d, input$EDbrush)
        return(d)
    })

    output$combEffect <- renderPlot({
        d <- dataEDplot()
        meth <- input$combMeth
        g <- ggplot(d, aes_string(x = "value", y = meth)) + geom_point()
        s <- input$dataCombEffect_rows_selected
        if (length(s)) {
            d2 <- dataEDplotBrushed()
            g <- g + geom_point(aes_string(x = "value", y = meth), data = d2[s, , drop = FALSE], color = "red", size = 5)
        }
        g <- g + xlab("Response") + ylab(input$combMeth)
        g <- g + theme(plot.background = element_blank(),legend.background = element_blank())
        return(g)
    }, bg = "transparent")

    getCITable <- reactive({
        d <- dataEDplotBrushed()
        d <- d[, c(drugA(), drugB(), "value", "HSA", "Bliss", "LoeweExcess", "Chou")]
        d <- round(d, digits = 3)
        colnames(d) <- c(drugA(), drugB(), "Response", "HSA", "Bliss", "LoeweExcess", "Chou")
        return(d)
    })

    output$dataCombEffect <- DT::renderDataTable({
        d <- getCITable()
        return(DT::datatable(d, rownames = TRUE,  options = list(dom = 'lrtip', autowidth = FALSE)))
    })

    output$tableDL <- downloadHandler(filename=function(){paste0(drugA(), "_", drugB(), ".csv")}, content = function(file){
        if (file.exists(file)) file.remove(file)
        d <- getCITable()
        write.csv(d, file, row.names = FALSE, quote = FALSE)
    })


    ## exported image with results
    output$isobolUI <- renderUI({
        if ("isobol" %in% input$reportChoices) {
        ui <- fluidRow(column(width = 4, selectInput("isobolTypeReport", "Isobologram type", choices = c("Linear" = "linear", "Steel and Peckham" = "SP"))), column(width = 4, numericInput("effectIsoReport", "Effect", value = 0.5, min = 0, max = 1, step = 0.1)))
        } else ui <- NULL
        return(ui)
    })

    reportPlot <- reactive({
        object <- object()
        choices <- input$reportChoices
        glist <- list()
        if ("Data" %in% choices) {
            g1 <- plot(object, type = "heatmap" , what = "value") + ggtitle("Data")
            glist <- c(glist, g1 = list(g1))
        }
        if ("HSA" %in% choices) {
            g2 <- plot(object, type = "heatmap" , what = "HSA") + ggtitle("HSA")
            glist <- c(glist, list(g2))
        }
        if ("Bliss" %in% choices)  {
            g3 <- plot(object, type = "heatmap" , what = "Bliss") + ggtitle("Bliss")
            glist <- c(glist, list(g3))
        }
        if ("indResp" %in% choices) {
            gA <- plot(object, type = "ind", drug = drugNames(object)[1])
            gB <- plot(object, type = "ind", drug = drugNames(object)[2])
            glist <- c(glist, list(cowplot::plot_grid(gA, gB, ncol = 1)))
        }
        if ("isobol" %in% choices) {
            g5 <- isobologram(object, effect = input$effectIsoReport, mode = input$isobolTypeReport) + ggtitle(paste("Isobologram, effect", input$effectIsoReport))
            glist <- c(glist, list(g5))
        }
        if ("LoewePred" %in% choices) {
            g6 <- plot(object, type = "heatmap", what = "Loewe") + ggtitle("Loewe predicted")
            glist <- c(glist, list(g6))
        }
        if ("LoeweExcess" %in% choices) {
            g7 <- plot(object, type = "heatmap", what = "LoeweExcess") + ggtitle("Loewe excess")
            glist <- c(glist, list(g7))
        }
        if ("Chou" %in% choices) {
            g8 <- plot(object, type = "heatmap", what = "Chou") + ggtitle("Chou")
            glist <- c(glist, list(g8))
        }
        p <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label(input$reportTitle, fontface = "bold", size = 20), cowplot::plot_grid(plotlist = glist, scale = 0.95), ncol = 1, rel_heights = c(0.05, 1))
        return(p)
    })

    output$downloadReport <- downloadHandler(filename = function(){
        paste(drugA(), drugB(), ".png", sep = "")
        },
        content = function(file){
            png(file, width = 1200, height = 800)
            print(reportPlot())
            dev.off()
        })

    output$reportTitleUI <- renderUI({
        object <- object()
        textInput("reportTitle", label = "Report title", value = paste(drugNames(object), collapse = " - "))
        })

    getPage <- function() {
        #        return(includeHTML("help.html"))
        return(includeMarkdown("help.Rmd"))
    }
    output$help <- renderUI({withMathJax(getPage())})
})




