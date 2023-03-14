#
# Código inicial, modelo matemático e idea:
# Michael Höhle, 'Flatten the COVID-19 curve', Mar 16, 2020
# https://staff.math.su.se/hoehle/blog/2020/03/16/flatteningthecurve.html
# Desarrollado por Tinu Schneider tinu@tinuschneider.ch 
# Code is on Github github.com/tinu-schneider/Flatten_the_Curve 





library(shiny)
library(gridExtra)
library(tidyverse)
library(deSolve)
library(lubridate)
library(stats)


## Configuration
theme_set(theme_minimal(base_size = 18))


## CONSTANTS

# Population size 
N <- 654214 

# Rate at which person stays in the infectious compartment (disease specific and tracing specific)
gamma <- 1/12

# Initial number of infected people

casos_al_inicio_confinamiento<-read.csv("https://raw.githubusercontent.com/ibaitamayo/infectious/main/COVID_19/Covid_data/ccaa_covid19_confirmados_pcr_long.csv",encoding = "UTF-8")%>%filter(CCAA=="Navarra" & fecha=="2020-03-20")%>%dplyr::select(total)


# casos_al_inicio_confinamiento<-read.csv("https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/ccaa_covid19_confirmados_pcr_long.csv",encoding = "UTF-8")%>%filter(CCAA=="Navarra" & fecha=="2020-03-20")%>%select(total)
# se seleccionan los positivos de 5 días tras confinamiento (mediana hasta  ser clasificado como positivo), ya que los positivos de un momento se infectaron 5 días atrás

#https://www.mscbs.gob.es/ciudadanos/ene-covid/docs/ESTUDIO_ENE-COVID19_PRIMERA_RONDA_INFORME_PRELIMINAR.pdf 
seroprevalencia<-0.058

Infectados_totales_en_fecha_Seroprevalencia<-round(0.058*N,0)


# casos_al_medir_seroprevalencia<-read.csv("https://raw.githubusercontent.com/datadista/datasets/master/COVID%2019/ccaa_covid19_confirmados_pcr_long.csv",encoding = "UTF-8")%>%filter(CCAA=="Navarra" & fecha=="2020-05-01")%>%select(total)

casos_al_medir_seroprevalencia<-read.csv("https://raw.githubusercontent.com/ibaitamayo/infectious/main/COVID_19/Covid_data/ccaa_covid19_confirmados_pcr_long.csv",encoding = "UTF-8")%>%filter(CCAA=="Navarra" & fecha=="2020-05-01")%>%dplyr::select(total)
#se seleccionan 10 días atrás, (tiempo medio para seroconversión desde la infección)
porcent_positivos_testados<-casos_al_medir_seroprevalencia/Infectados_totales_en_fecha_Seroprevalencia

n_init <- as.numeric(round(casos_al_inicio_confinamiento /porcent_positivos_testados))
# n_init <- 1830

# ccaalist<-as.list(levels(as.factor(x$CCAA)))
# colnames(x)[grep(pattern="fecha",x=colnames(x))]<-"fecha"

# # Grid where to evaluate
# max_time <- 365 # 150
# times <- seq(0, max_time, by = 0.1)

# R0 for the beta and gamma values
# R0 <- beta*N/gamma

# calculate beta
# Infectious contact rate - beta = R0/N*gamma and when R0  ~2.25 then  2.25/N*gamma
# beta <- 4.5e-07 





# Function to compute the derivative of the ODE system
# -----------------------------------------------------------
#  t - time
#  y - current state vector of the ODE at time t
#  parms - Parameter vector used by the ODE system (beta, gamma)

sir <- function(t, y, parms, 
                social_dist_period, 
                reduction,delay=c(12,12,2,5)) {
    
    beta0 <- parms[1]
    gamma <- parms[2]
    
    # Reduce contact rate 
    beta_z = if_else(t <= social_dist_period[1], beta0, if_else(t <= as.integer(social_dist_period[1]+delay[1])& t > as.integer(social_dist_period[1]),beta0*1+(as.integer(t-social_dist_period[1])*((beta0*reduction[1]-beta0)/delay[1])), if_else(t <= social_dist_period[2] & t>as.integer(social_dist_period[1]+delay[1]), beta0 * reduction[1],if_else(t>social_dist_period[2]& t <= as.integer(social_dist_period[2]+delay[2]), beta0*reduction[1]+(as.integer(t-social_dist_period[2])*((beta0*reduction[2]-beta0*reduction[1])/delay[2])),beta0 * reduction[2]))))#modifico la función para adaptación progresiva  
    beta_t = if_else(t <= social_dist_period[3], beta_z, if_else(t <= as.integer(social_dist_period[3]+delay[3])& t > as.integer(social_dist_period[3]),beta_z*1+(as.integer(t-social_dist_period[3])*((beta0-(beta0*reduction[3]))/delay[3])), if_else(t <= social_dist_period[4] & t>as.integer(social_dist_period[3]+delay[3]), beta0 * reduction[3],if_else(t <= as.integer(social_dist_period[4]+delay[4]), beta0*reduction[3]+(as.integer(t-social_dist_period[4])*((beta0*reduction[4]-beta0*reduction[3])/delay[4])),beta0 * reduction[4]))))#modifico la función para adaptación progresiva  
    S <- y[1]
    I <- y[2]
    
    return(list(c(S = -beta_t * S * I, 
                  I =  beta_t * S * I - gamma * I)))
}


## we assume that some globals exist...
solve_ode <- function(sdp, red, typ, beta, max_time) {
    
    # Grid where to evaluate
    times <- reactive({ seq(0, max_time, by = 0.1) })
    
    ode_solution <- lsoda(y = c(N - n_init, n_init), 
                          times = times(), 
                          func  = sir, 
                          parms = c(beta, gamma), 
                          social_dist_period = sdp,
                          reduction = red) %>%
        as.data.frame() %>%
        setNames(c("t", "S", "I")) %>%
        mutate(beta = beta, 
               gama = gamma,
               R0 = N * beta / gamma, 
               s  = S / N, 
               i  = I / N, 
               type = typ)
    
    daily <- ode_solution %>%
        filter(t %in% seq(0, max_time, by = 1)) %>%
        mutate(C = if_else(row_number() == 1, 0, lag(S) - S), 
               c = C / N)
    
    daily
}

solve_ode2 <- function(sdp, red, typ, beta, diashastaconfi) {
    
    # Grid where to evaluate
    times <- reactive({ seq(0, diashastaconfi, by = 0.1) })
    
    ode_solution <- lsoda(y = c(N - 1, 1), 
                          times = times(), 
                          func  = sir, 
                          parms = c(beta, gamma), 
                          social_dist_period = diashastaconfi,
                          reduction = c(0,0)) %>%
        as.data.frame() %>%
        setNames(c("t", "S", "I")) %>%
        mutate(beta = beta, 
               gama = gamma,
               R0 = N * beta / gamma, 
               s  = S / N, 
               i  = I / N)
    
    daily <- ode_solution %>%
        filter(t %in% seq(0, diashastaconfi, by = 1)) %>%
        mutate(C = if_else(row_number() == 1, 0, lag(S) - S), 
               c = C / N)
    
}


discreto<- function(time_ini=0,time_end=-300,beta=beta,reduction=c(1,1),social_dist_period=0,delay=c(0,0),I_ini=n_init,S_ini=N,R_ini=0){
    tvec <- seq(time_ini,time_end,-1) 
    I_ini->I
    S_ini->S
    R_ini->R
    beta->beta0
    beta_t<-beta0
    fecha2<-NA
    ts = data.frame(cbind(tvec,S,I,R,beta0,fecha2))
    colnames(ts) = c("time","S","I","R","beta","fecha2") 
    ct=1 #a counter to index array 
    for (t in tvec) 
    {
        ts[ct,] = c(t,S,I,R,beta0,fecha2) 
        Sp = S - 1*(-beta_t*S*I)
        Ip = I - 1*(+beta_t*S*I -gamma*I)
        Rp = R - 1*(+gamma*I)
        fecha2=as.character(ymd("2020-03-15")+lubridate::days(t))
        S = Sp 
        I = Ip 
        R = Rp 
        ct = ct + 1 
    } #finish loop 
    return(ts)
}







# add results with intervention and plot
run <- function(sdp, red, r0, max_time) {
    
    beta <- r0 / N * gamma
    
    ode_solution_daily <- solve_ode(
        sdp = c(0, max_time),  # social_dist_period
        red = c(1, 1),         # reduction
        typ = "sin", 
        beta = beta, 
        max_time
    )
    
    # solve with interventions
    ode_solution2_daily <- solve_ode(
        sdp = sdp,
        red = red,
        typ = "con", 
        beta = beta, 
        max_time
    )
    
    # Combine the two solutions into one dataset
    ode_df <- rbind(ode_solution_daily, ode_solution2_daily)%>%mutate(fecha=ymd("2020-03-15")+days(t))
    
    
}
run2 <- function(r0) {
    beta <- r0 / N * gamma
    discreto(beta = beta)->fecha2
    
    return(fecha2)
}

run3<-function(r0){
    beta <- r0 / N * gamma
    
}

plot_result <- function(ode_df, sdp, max_time, y_axis_fixed,fecha2) {    
    
    # The final size in the two cases:
    final_sizes <- ode_df %>%
        group_by(type) %>%
        filter(row_number() == n()) %>%
        mutate("fracción final" = scales::percent(1 - s, accuracy = 1)) %>%
        dplyr::select("fracción final", intervenciones = type) %>% 
        arrange(desc(intervenciones))
    
    # Plot
    if (y_axis_fixed) {
        y_max <- 0.09
    } else {
        y_max <- max(ode_df$C[ode_df$type=="con"], na.rm = TRUE) * 1.05
    }
    y_arrow <- y_max * 0.975
    y_text  <- y_arrow + y_max * 0.01 
    col_sdp <- "lightblue3"
    satura<-max(ode_df$C[ode_df$type=="con" &ode_df$t>sdp[1]& ode_df$t<sdp[2]])
    diasatura<-as.character(ode_df%>%filter(type=="con" & t>sdp[2]& c>satura)%>%filter(t==min(t))%>%dplyr::select(t))
    fechasatura<-ifelse(diasatura=="numeric(0)","",as.character(ymd("2020-03-15")+days(as.numeric(diasatura))))
    
    x_labs <- sort(c(0, 100, 200, 300, 365, sdp))
    
    
    
    pp <- ggplot(ode_df, 
                 aes(x = t, 
                     y = 0, 
                     xend = t, 
                     yend = C, 
                     color = type)) + 
        geom_segment(data=ode_df%>%filter(type=="con"),alpha = 0.7) + 
        geom_line(aes(x = t, y = C)) + 
        geom_point(aes(x = t, y = C,colour=type)) + 
        labs(
            x = "Días", 
            y = NULL, 
            subtitle = "Nuevos casos diarios", 
            caption  = paste("Fecha estimada de primera infección : \n",as.character(fecha2%>%filter(I<=1 & I>0)%>%filter(I==max(I))%>%dplyr::select(fecha2)),sep="")) +
        # facet_grid(scales = "free")+
        # scale_x_discrete(labels = ode_df$fecha[1:t,], 
        #                    breaks = x_labs) +
        scale_x_continuous(labels = as.character(ymd("2020-03-15")+days(ode_df$t[ode_df$type=="con"]))[seq(0,y_max,by=30)],breaks = ode_df$t[ode_df$type=="con"][seq(0,y_max,by=30)])+
        # scale_x_date(labels = as.Date(ode_df$fecha[ode_df$type=="con"]),breaks = as.Date(ode_df$fecha[ode_df$type=="con"]))+
        # scale_y_continuous(labels = scales::percent,
        #                    limits = c(0, y_max)) +
        theme(panel.grid.minor = element_blank()) +
        scale_color_brewer(name = "Intervenciones", 
                           type = "qual", 
                           palette = 1, 
                           guide = guide_legend(reverse = TRUE)) +
        # sdp 
        geom_vline(xintercept = c(sdp), lty = 2, color = "darkgray") +
        # geom_vline(xintercept = c(sdp), lty = 3, color = "darkgray") +
        geom_hline(yintercept=satura)+
        geom_hline(yintercept=ode_df$C[ode_df$type=="con" &ode_df$t==sdp[3]],lty=2)+
        # geom_vline(xintercept = ode%>%filter(type=="con")%>%mutate(dato=if_else(t<sdp[2])),
        #                min(ode_df$t[ode_df$C>max(ode_df$C[ode_df$type=="con" & ode_df$t<sdp[2]])&ode_df$t[ode_df$type=="con"]>sdp[2]]), lty = 2, color = "darkgray") +
        geom_text(aes(x = sdp[1] + (sdp[2] - sdp[1])/2, 
                      y = y_text, 
                      label = "Confinamiento"),
                  vjust = 0,
                  color = col_sdp) +   
        geom_text(aes(x = as.numeric(diasatura)-15,
                      y = satura * 1.05,
                      label =fechasatura),
                  vjust = 0,
                  color = "grey") +
        geom_text(aes(x = sdp[3] + (sdp[4] - sdp[3])/2, 
                      y = y_text, 
                      label = "Descontrol"),
                  vjust = 0,
                  color = col_sdp) +   
        geom_segment(aes(
            x = sdp[1],
            y = y_arrow, 
            xend = sdp[2] * 0.99, # shorten
            yend = y_arrow
        ),
        size = 0.3, 
        color = col_sdp,
        arrow = arrow(length = unit(2, "mm"))) +
        geom_segment(aes(
            x = sdp[3], # shorten
            y = y_arrow, 
            xend = sdp[4],
            yend = y_arrow
        ),
        size = 0.3, 
        color = col_sdp,
        arrow = arrow(length = unit(2, "mm")))  +
        # Add final size as table
        # annotation_custom(tableGrob(final_sizes, 
        #                             rows = NULL,
        #                             theme = ttheme_minimal(
        #                                 core    = list(fg_params = list(hjust = 0, x = 0.1)),
        #                                 rowhead = list(fg_params = list(hjust = 0, x = 0))
        #                             )),
        #                   xmin = sdp[2],
        #                   xmax = sdp[2]*1.2,
        #                   ymin = y_max * 0.5,
        #                   ymax = y_max * 0.6
        # )+
        theme(legend.position="none")
    
    print(pp)
}


ui <- fluidPage(
    
    titlePanel("Modelizando COVID19 e intervenciones no farmacológicas en Navarra"),
    
    sidebarLayout(
        sidebarPanel(
            sliderInput("r0", 
                        div(HTML("Valor de R<sub>0</sub>")), 
                        min   = 0.5, 
                        max   = 15.0, 
                        value = 2.8, 
                        step  = 0.1),
            checkboxInput("y_axis_fixed",
                          div(HTML("Fija eje y  (útil para comparar diferentes R<sub>0</sub>)")),
                          value = FALSE),
            br(),
            sliderInput("x_max", 
                        "Max de dias para el modelo", 
                        min   = 45, 
                        max   = 365, 
                        value = 250,
                        step  =  15),
            hr(),
            br(),
            sliderInput("sdp",
                        "Confinamiento (pds 1)",
                        min   =   0,
                        max   = 150,
                        value = 50, 
                        step  =   5), 
            
            sliderInput("red_one",
                        "Reducción contactos confinamiento (pds 1)",
                        min   = 0.0,
                        max   = 1.0,
                        value = 0.2, 
                        step  = 0.05), 
            sliderInput("red_two",
                        "Reducción contactos posterior (pds 2)",
                        min   = 0.0,
                        max   = 1.0,
                        value = 0.35, 
                        step  = 0.05), 
            br(), br(),
            sliderInput("sdp_unc",
                        "Uncontrolled period ",
                        min   =   0,
                        max   = 150,
                        value = c(113,121), 
                        step  =   5), 
            br(), br(),
            sliderInput("red_three",
                        "Reducción de contactos en periodo descontrolado (Uncontrolled)",
                        min   = 0.0,
                        max   = 1.0,
                        value = 0.65, 
                        step  = 0.05), 
            sliderInput("red_four",
                        "Reducción contactos tras periodo descontrolado (Gaining control)",
                        min   = 0.0,
                        max   = 1.0,
                        value = 0.35, 
                        step  = 0.05), 
            
            br(),
            
            
            img(src="data:image/png;base64,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"),
            hr(),
            tags$div("Unidad de Metodología, 2020"), 
            tags$a(href="http://www.navarrabiomed.es/es/servicios/unidades-de-apoyo/metodologia"),
            br(),
            br(),
            
            tags$div("Adaptado del codigo en Github:"),
            # br(),
            tags$a(href="https://github.com/tinu-schneider/Flatten_the_Curve", 
                   "github.com/tinu-schneider/Flatten_the_Curve"),
            br(), br(),
            tags$div("Trabajo bajo licencia "),
            tags$a(href="http://creativecommons.org/licenses/by-sa/4.0/", 
                   "Creative Commons Attribution-ShareAlike 4.0 International License") 
            
        ),
        
        mainPanel(
            # p("Note: This tool is not intended to create a prediction."),
            plotOutput("chart", height = "500px"), 
            br(),
            hr(),
            h4("...a hombros de gigantes:" ),
            tags$a(href="https://staff.math.su.se/hoehle/blog/2020/03/16/flatteningthecurve.html", 
                   "Michael Höhle y Tinu Schneider, 'Flatten the COVID-19 curve'"), 
            p(),
            br()
        )
    )
)


server <- function(input, output) {
    
    
    res <- reactive({
        run(sdp = c(0,input$sdp,input$sdp_unc), 
            red = c(input$red_one, input$red_two,input$red_three, input$red_four), 
            r0  = input$r0, 
            max_time = input$x_max)
    })
    res2<-reactive({
        run2(r0  = input$r0)
    })
    
    # print(res())
    output$chart <- renderPlot({
        plot_result(res(), c(0,input$sdp,input$sdp_unc), input$x_max, input$y_axis_fixed,res2())
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
