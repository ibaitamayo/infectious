
library(leaflet)    
library(shiny)
library(shinydashboard)
library(rgdal)
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(lubridate)



# y<-read.csv("https://raw.githubusercontent.com/ibaitamayo/infectious/main/COVID_19/Covid_data/ccaa_covid19_test_realizados.csv",encoding = "UTF-8")
y<-read.csv("https://raw.githubusercontent.com/ibaitamayo/infectious/main/COVID_19/Covid_data/ccaa_vacunas.csv",encoding = "UTF-8")
y<-y[,-c(1,4:9,11:12)]
colnames(y)[grep(pattern="fecha",x=tolower(colnames(y)))]<-"fecha"
codcom<-y%>%distinct(CCAA,cod_ine)

# seroprev1_2<-read.csv("https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/provincias_estudio_prevalencia_anticuerpos_primera_y_segunda_ronda.csv",encoding = "UTF-8",dec = ",")
# 
# pob2019 <- read_delim("C:/Users/USER-PC/Downloads/pob2019.csv", ";", escape_double = FALSE, trim_ws = TRUE)[-1,]%>%mutate(CP=as.integer(strtrim(Provincias,width = 2)),poblacion=gsub(Total,pattern = "\\.",replacement = ""))%>%select(CP,poblacion)
# prov_cod_km2 <- read_excel("C:/Users/USER-PC/Downloads/prov_cod_km2.xls", 
#                            sheet = "list-pro", col_types = c("text", 
#                                                              "skip", "numeric", "numeric", "skip", 
#                                                              "skip", "skip"))
# pob2019%>%full_join(prov_cod_km2%>%mutate(CP=as.integer(CP)),by="CP")%>%mutate(densidad=as.numeric(poblacion)/Superficie)->Densidad
# 
# 
# Densidad%>%full_join(seroprev1_2%>%mutate(CP=1:n())%>%select(CP,total_primera_ronda,total_porcentaje_primera_ronda,porcentaje_seroconversion_IgG.,total_porcentaje_segunda_ronda),by="CP")->serodensidad
# 
# x<-read.csv("https://raw.githubusercontent.com/ibaitamayo/infectious/main/COVID_19/Covid_data/ccaa_covid19_confirmados_pcr_long.csv",encoding = "UTF-8")
# gg<- read.csv("https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/ccaa_covid19_datos_isciii_nueva_serie.csv",encoding = "UTF-8")


gg<- read.csv("https://raw.githubusercontent.com/ibaitamayo/infectious/main/COVID_19/Covid_data/ccaa_covid19_datos_isciii_nueva_serie.csv",encoding = "UTF-8")


colnames(gg)[grep(pattern="fecha",x=colnames(gg))]<-"fecha"
colnames(gg)[grep(pattern="ccaa",x=colnames(gg))]<-"CCAA"
ccaalist<-as.list(levels(as.factor(gg$CCAA)))
gg%>%group_by(CCAA)%>%arrange(fecha)%>%mutate(diarios=num_casos_prueba_pcr)%>%arrange(fecha)%>%mutate(num_dias=1:n())->x
# unique(x$num_dias[x$fecha=="2020-03-15"])->confi
x$fecha2<-ymd(x$fecha)

# View(read.csv("https://raw.githubusercontent.com/datadista/datasets/master/COVID 19/provincias_estudio_prevalencia_anticuerpos_primera_y_segunda_ronda.csv",encoding = "UTF-8"))

movingAverage <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}


# x<-read.csv("https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/ccaa_covid19_confirmados_pcr_long.csv",encoding = "UTF-8")
# ccaalist<-as.list(levels(as.factor(x$CCAA)))
# colnames(x)[grep(pattern="fecha",x=colnames(x))]<-"fecha"
# x%>%group_by(CCAA)%>%arrange(fecha)%>%mutate(diarios=total-shift(total,n=1,type="lag",fill=0))%>%arrange(fecha)%>%mutate(num_dias=1:n())->x

codcom<-x%>%distinct(CCAA,cod_ine)
y<-left_join(codcom,y%>%dplyr::select(-CCAA),by="cod_ine")

y%>%group_by(CCAA)%>%arrange(fecha)%>%mutate(porcent_dosis_por_100_habit=as.numeric(gsub(pattern=",",replacement = ".",x = Porcentaje.de.dosis.administradas.por.100.habitantes)),POrcent_pauta_completa=as.numeric(gsub(pattern=",",replacement = ".",x = Porcentaje.con.pauta.completa)),fecha2=dmy(fecha))->y
Codigo<-c("ES.AN","ES.AR","ES.AS","ES.IB","ES.CN","ES.CB","ES.CM","ES.CL","ES.CT","ES.CE","ES.VC","ES.EX","ES.GA","ES.MD","ES.ML","ES.MC","ES.NC","ES.PV","ES.RI")
cod_ine<-c(1,2,3,4,5,6,8,7,9,18,10,11,12,13,19,14,15,16,17)
Codigo1<-cbind(data.frame(cod_ine),data.frame(Codigo))
codif<-full_join(codcom,as.data.frame(Codigo1),by="cod_ine")


mapas_comunidades<-readOGR("https://raw.githubusercontent.com/ibaitamayo/utilidades/main/Mapas/spain-comunidad-with-canary-islands.json")
mapas_comunidades@data%>%mutate(HASC_1=as.character(HASC_1))%>%
  mutate(HASC_1=ifelse(HASC_1=="ES.NA","ES.NC",ifelse(HASC_1=="ES.PM","ES.IB",ifelse(HASC_1=="ES.LO","ES.RI",ifelse(HASC_1=="ES.MU","ES.MC",HASC_1)))))%>%
  mutate(Codigo=HASC_1)%>%left_join(codif,by="Codigo")%>%
  dplyr::select(Codigo,CCAA,cod_ine)->mapas_comunidades@data
mapas_comunidades@data$id<-as.factor(mapas_comunidades@data$cod_ine)


ui <- fluidPage(verticalLayout(
    # selectInput(inputId = "CCAA", label="Elige la comunidad:",
    #             choices=ccaalist),
    
    # mainPanel(
        
    # place the contents inside a box
    shinydashboard::box(
        width = 12
        , title = "Selecciona la Comunidad Autónoma"
        # separate the box by a column
        , column(
            width = 4
            , shiny::actionButton( inputId = "clearHighlight"
                                   , icon = icon( name = "eraser")
                                   , label = "Borrar selección"
                                   , style = "color: #fff; background-color: #D75453; border-color: #C73232"
            )
        )
        # separate the box by a column
        , leaflet::leafletOutput( outputId = "myMap"
                                      
            
        )
    ), # end of the box
    plotOutput("distPlot")
   
    )
) # end of fluid page

# create the server
server <- function( input, output, session ){

    # function to create foundational map
    foundational.map <- function(){
        leaflet() %>%
            setView( lng = -3
                    , lat = 40
                    , zoom = 5 ) %>%
            addProviderTiles("Esri.WorldShadedRelief") %>%
            addPolygons( data = mapas_comunidades
                         , fillOpacity = 0
                         , opacity = 0.2
                         , color = "#000000"
                         , weight = 2
                         , layerId = mapas_comunidades$CCAA
                         , group = "click.list") 
    }
    
    # reactiveVal for the map object, and corresponding output object.
    myMap_reval <- reactiveVal(foundational.map())
    
    output$myMap <- renderLeaflet({
        myMap_reval()
    }) 
    
    # To hold the selected map region id.
    click.list <- shiny::reactiveValues( ids = vector() )
    
    shiny::observeEvent( input$myMap_shape_click, ignoreNULL = T,ignoreInit = T, {

        # If already selected, first remove previous selection
        if(length(click.list)>0)
        {
            remove_id = click.list$ids
            lines.of.interest <- mapas_comunidades[ which( mapas_comunidades$CCAA %in% remove_id) , ]
            leaflet::leafletProxy( mapId = "myMap" ) %>%
                addPolylines( data = lines.of.interest
                              , layerId = lines.of.interest@data$id
                              , color = "#000000"
                              , weight = 2
                              , opacity = 0.2) 

        }
        
        # add current selection
        click <- input$myMap_shape_click
        click.list$ids <- click$id  # we only store the last click now!
        lines.of.interest <- mapas_comunidades[ which( mapas_comunidades$CCAA %in% click.list$ids ) , ]
        
        # print(click)

        if( is.null( click$id ) ){
            req( click$id )
        } else if( !click$id %in% lines.of.interest@data$id ){
            leaflet::leafletProxy( mapId = "myMap" ) %>%
                addPolylines( data = lines.of.interest
                              , layerId = lines.of.interest@data$id
                              , color = "#6cb5bc"
                              , weight = 5
                              , opacity = 1
                ) 
  
        }

    }) # end of shiny::observeEvent({})
    
    # oberver for the clearHighlight button.
    shiny::observeEvent( input$clearHighlight, {
        click.list$ids <- NULL
        myMap_reval(foundational.map()) # reset map.
    }) 
    output$distPlot <- renderPlot({
      
      if (length(click.list$ids)>0){

        filtered<-x%>%filter(CCAA==click.list$ids)%>%filter(diarios>=0)
        pt2<-ggplot(data=filtered,aes(x=fecha2,y=diarios))+geom_point()+ggtitle(paste("Casos diarios en",click.list$ids,sep=" "))+geom_line(aes(y=movingAverage(x = diarios,n = 7,centered = TRUE),colour="red"))+ylab("Casos diarios")+xlab("fecha")+geom_vline(xintercept = ymd("2020-03-15"))+theme_light()+theme(legend.position="none")
        filtered2<-y%>%filter(CCAA==click.list$ids)
        pt1<-ggplot(data=filtered2,aes(x=fecha2))+geom_line(aes(y=POrcent_pauta_completa), colour = "blue")+ggtitle(paste("Porcentaje de vacunas en",click.list$ids,sep=" "))+ylab("pauta completa (%)")+xlab("fecha")+ylim(0,100)+theme(legend.position="none")+theme_light()#+geom_line(aes(y=porcent_dosis_por_100_habit,color="Porcentaje de dosis por cada 100 hab"))
        grid.arrange(pt1,pt2)
      }
      

    })
    
} 

shiny::shinyApp( ui = ui, server = server)

