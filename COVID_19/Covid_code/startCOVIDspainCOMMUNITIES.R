#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(data.table)
x<-read.csv("https://raw.githubusercontent.com/ibaitamayo/infectious/main/COVID_19/Covid_data/ccaa_covid19_confirmados_pcr_long.csv",encoding = "UTF-8")
ccaalist<-as.list(levels(as.factor(x$CCAA)))
colnames(x)[grep(pattern="fecha",x=colnames(x))]<-"fecha"
# Define UI for application that draws a plot
ui <- fluidPage(
    
    selectInput(inputId = "CCAA", label="Elige la comunidad:",
                choices=ccaalist),
    mainPanel(
        plotOutput("distPlot")
    )
)


# Define server logic required to draw
server <- function(input, output) {


    output$distPlot <- renderPlot({
        filtered<-x%>%group_by(CCAA)%>%arrange(fecha)%>%mutate(diarios=total-data.table::shift(total,n=1,type="lag",fill=0))%>%arrange(fecha)%>%mutate(num_dias=1:n())%>%filter(CCAA==input$CCAA)%>%filter(diarios>=0)

        ggplot(data=filtered,aes(x=num_dias,y=diarios))+geom_point()+ggtitle("PCRs Positivas por comunidad autónoma")+geom_smooth(method="gam")+geom_vline(xintercept = c(24,74,80,95))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)



