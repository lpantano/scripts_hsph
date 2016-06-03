shinyUI(pageWithSidebar(
    headerPanel("Plot Expression"),
    sidebarPanel(
        textInput("gene", "Gene ID:", "ENSMUSG00000000628"),
        textInput("mirna", "miRNA ID:", "mmu-let-7g-5p"),
        actionButton("do","Update View")   
    ),
    mainPanel(
        plotOutput("distPlot")
    )
))