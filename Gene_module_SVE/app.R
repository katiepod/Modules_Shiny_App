ui <- fluidPage(
  h1("Calculation of module activity"),
  h3(tags$div(
    "Follow this link to learn more about each module:", 
    tags$br(),
    a("Module Annotation Pages",
      href = "http://10.11.14.33/modules/m1.html"))),
  p(style = "font-family:Times;color:gray;size=12", "If you need to convert gene or protein IDs into Affymetrix Probeset IDs, you can do so here: ",
    a("DAVID",
      href = "https://david.ncifcrf.gov/conversion.jsp")),
  fluidRow(
    column(3,
           fileInput("file1", 
                     "Upload a text (.txt) file \n (Affymetrix probeset IDs in the first column, fold changes in the second column)",
                     accept = c(
                       "text/csv",
                       "text/comma-separated-values,text/plain",
                       ".csv")),
           textInput(inputId = "title", 
                     label = "Give a name to your experiment",
                     value = "Example: daf-2 vs. WT"),
           actionButton(inputId = "go", label = "Go")
    )
  ),
  
  fluidRow(
    column(3,
           tableOutput("example")),
    column(4,
           tableOutput("Projection")),
    column(5,
           plotOutput("projectionPlot"))
  )
)

server <- function(input, output) {
  
  output$example <- renderTable({
    daf2 <- read.table("daf-2_FCs.txt", sep="\t")
    mm1 <- as.matrix(daf2[,2])
    rownames(mm1) <- daf2[,1]
    colnames(mm1) <- "  "
    head(mm1)
  }, caption="What your input file should look like:",caption.placement = getOption("xtable.caption.placement", "top")
  )
  
  
  SP <- eventReactive(input$go, {
    inFile <- input$file1
    X <- as.matrix(read.table(inFile$datapath, sep = "\t",row.names = 1))
    X <- (X-mean(X))/sd(X)
    
    # Upload module definition matrix:
    load("GenModules.Affy.S.matrix.Rdata")
    
    # Define scalar projection function:
    getScalarProjection = function(X, S) {
      # Make sure both are matrices and keep probesets that are present in both
      stopifnot(class(X) == "matrix", class(S) == "matrix")
      commonRows <- intersect(rownames(X), rownames(S))
      X <- X[commonRows, ]
      S <- S[commonRows, ]
      
      # Name S columns if not named already
      if(is.null(colnames(S))) colnames(S) = as.character(1:ncol(S))
      
      # Get dot product
      d.p = t(S) %*% X
      
      # Find magnitude (L2 norm) of each column of S
      S.norms = apply(S, 2, vector.norm)
      
      # Divide dot product by norms
      A = d.p / S.norms
      colnames(A) = colnames(X)
      rownames(A) = colnames(S)
      
      return(A)
    }  
    vector.norm = function(x) {sqrt(sum(x^2))}
    
    # Define Z score function:
    getProjectionZscores = function(X, S, n = 100) {
      # Returns a list with z-scores for each
      # module in S for the given set of fold changes, X
      
      # X should be a one column matrix
      stopifnot(dim(X)[2] == 1)
      
      # Get true projection
      actual.proj = getScalarProjection(X, S)
      
      # Create random projections
      starting.matrix = matrix(X[,1], nrow = nrow(X), ncol = n)
      rand.matrix = apply(starting.matrix, 2, FUN = function(x) {
        sample(x, size = length(x), replace = FALSE)
      })
      rownames(rand.matrix) = rownames(X)
      rand.proj = getScalarProjection(rand.matrix, S)
      
      # Get distribution parameters for random projections
      rand.means = apply(rand.proj, 1, mean)
      rand.sd = apply(rand.proj, 1, sd)
      z.scores = (actual.proj[,1] - rand.means) / rand.sd
      
      # Set names
      names(z.scores) = colnames(S)
      
      return(z.scores)
    }
    
    # Define p-value function:
    convertZtoPvalues = function(z) {
      # Converts z-scores to p-values
      # Returned values are -log10 of p-values
      # Signs are preserved, e.g. negative z scores are given a negative p-values
      too.high = abs(z) > 37
      if(sum(too.high) > 0) {
        message("Cannot estimate p-values for extreme Z-scores (abs(z) > 37); setting to +/- 37...")
        z[z < -37] = -37
        z[z > 37] = 37
      }
      pvalue2sided=2*pnorm(-abs(z))
      log.p = -log(pvalue2sided, base = 10)
      find.negs = z < 0
      log.p[find.negs] = log.p[find.negs] * -1
      
      return(log.p)
    }
    
    # Define SVE function:
    getSignedVarianceExplained = function(X, S) {
      # Get scalar projection
      A = getScalarProjection(X, S)
      # Preserve signs
      signs = ((A[, 1] < 0) * -2) + 1
      dataPower <- numeric(nrow(A))
      for (i in 1:nrow(A)) dataPower[i] <- sum(A[i, ]^2)
      dataPower <- dataPower/sum(dataPower)
      names(dataPower) <- 1:nrow(A)
      dataPower = dataPower * signs
      return(dataPower)
    }
    
    # Calculate z-scores:
    z <- getProjectionZscores(X,S, n=100)
    
    # Calculate p-values:
    p <- convertZtoPvalues(z)
    p <- cbind (ModuleNum = as.numeric(colnames(S)), p.value = p)
    p <- data.frame(p)
    p$ModuleNum <- rownames(p)
    colnames(p) <- c("ModuleNum","-log10(p-value)")
    # Remove (-) sign from p-values
    p[,2] <- abs(p[,2])
    
    # Calculate SVE and merge w/p-values:
    SVE <-getSignedVarianceExplained(X,S)
    SVE <- cbind (ModuleNum = as.numeric(colnames(S)), SVE = SVE)
    SVE <- data.frame(SVE)
    SVE <- merge(SVE,p,by="ModuleNum")
    SVE <- SVE[order(-SVE$`-log10(p-value)`),]
    # row.names(SVE) <- NULL
    
  })
  
  output$Projection <- renderTable({
    SVE <- SP()
    SVE$ModuleNum <- as.integer(SVE$ModuleNum)
    row.names(SVE) <- NULL
    SVE
  })
  
  output$projectionPlot <- renderPlot({
    
    SVE <- SP()
    sd3_above <- mean(SVE[,2]) + 3*(sd(SVE[,2]))
    sd3_below <- mean(SVE[,2]) - 3*(sd(SVE[,2]))
    
    sd2_above <- mean(SVE[,2]) + 2*(sd(SVE[,2]))
    sd2_below <- mean(SVE[,2]) - 2*(sd(SVE[,2]))
    
    signif <- SVE[abs(SVE[,3]) > 3,]
    
    lims <- max(abs(SVE[,2]))
    library(ggplot2)
    ggplot()+
      scale_x_continuous(limits=c(-5,225))+
      scale_y_continuous(limits=c(-lims,lims))+
      
      ggtitle(input$title)+
      
      geom_hline(aes(yintercept=sd3_above), colour="orange") +
      annotate("text", x = 0, y = sd3_above*1.2, label = "3 S.D.", color="orange",size=4)+
      geom_hline(aes(yintercept=sd3_below), colour="green")+
      annotate("text", x = 0, y = sd3_below*1.2, label = "3 S.D.", color="green",size=4)+
      geom_hline(aes(yintercept=sd2_above), colour="orange", lty="dashed") +
      annotate("text", x = 0, y = sd2_above*1.2, label = "2 S.D.", color="orange",size=4)+
      geom_hline(aes(yintercept=sd2_below), colour="green", lty="dashed")+
      annotate("text", x = 0, y = sd2_below*1.2, label = "2 S.D.", color="green",size=4)+
      
      geom_point(data=SVE, aes(x = ModuleNum,y = SVE), color="darkgray", size=3) + 
      geom_point(data=signif, aes(x = ModuleNum,y = SVE), color="black", size=4, alpha=0.7) + 
      geom_text(data=signif, aes(x = ModuleNum,y = SVE, label=ModuleNum), hjust=-0.5, vjust=-0.25, size=4) +         
      annotate("text", x = 10, y = lims*0.8, label = "p > 0.001", color="darkgray",size=5)+
      annotate("text", x = 10, y = lims*0.9, label = "p < 0.001", color="black",size=5)+
      
      xlab("Module Number") + 
      ylab("Signed Variance Explained\n")+
      
      theme_bw()+
      theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())   
  })
  
}

shinyApp(ui, server)