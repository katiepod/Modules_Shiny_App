ui <- fluidPage(
  h1("Welcome to the", em("C. elegans"),"gene module website!", style = "color:black;font-family:arial;text-align: center"),
  
  h3("If you would like to know which modules are active in your condition of interest and you have a list of gene fold changes (condition vs. control)", style = "color:gray;font-family:arial;text-align:center"),
  
  h3(a("enter here",
       href = "https://kpodshivalova.shinyapps.io/Gene_module_SVE/", style="color:lightgreen;font-weight:bold"), style = "font-family:arial;text-align:center"),
  
  br(),
  
  h3("If you have a gene of interest and would like to know which modules it belongs to, or if you have a list of several genes and would like to know if they are enriched in a particular module",style = "color:gray;font-family:arial;text-align:center"),
  
  h3(a("enter here",
       href= "https://kpodshivalova.shinyapps.io/Enrichment_analysis/", style="color:lightgreen;font-weight:bold"), style = "font-family:arial;text-align:center"),
  
  br(),
  
  h3(img(src="Old_worm.png", align="center", height = 200, width = 200),style = "text-align:center")
)

server <- function(input, output) {}

shinyAppDir(".")

shinyApp(ui = ui, server = server)