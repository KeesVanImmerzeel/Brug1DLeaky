library(shiny)
library(shinyBS)
library(pracma)
library(ggplot2)
library(dplyr)
library(plotly)

# Suppress dplyr::summarise info
options(dplyr.summarise.inform = FALSE)

#' Polder function, Chapter 3.1.2, formula 24, page 862.
#'
#' @param x First argument in the Polder function
#' @param y Second argument in the Polder function
#' @references Formula 24 page 862.
#' @return Result of the Polder function
#' @export
P <- function(x, y) {
      0.5*exp(2*x)*pracma::erfc(x/y+y) + 0.5*exp(-2*x)*pracma::erfc(x/y-y)
}

#' Calculate the head Phi in response to a 
#' Sudden drawdown of the surface water level, which is kept constant thereafter.
#' 
#' BI. One-dimensional groundwater flow
#' BI-2. One-dimensional groundwater flow in a semi-infinite field
#' The soil is assumed to be homogeneous
#' The flow is non-periodic
#' The boundary condition at x -- 0 a is given head or drawdown
#' Leaky aquifers with variable head at x = 0
#' Sudden drawdown of the surface water level, which is kept constant thereafter.
#' p. 66 Bruggeman Formula 123.32
Br_123_32 <- function(x, t, h, Eta, Labda) {
      Phi <- h * P(x / (2 * Labda), sqrt(Eta * t))
      return(Phi)
}

#' Calculate the Head change Phi [L] at distance x and time t due to a stress change h (= sudden drawdown) at t=0 at x=0
#'
#' @param x Distance [L]
#' @param t Time [T]
#' @param S Storage coefficient [-]
#' @param kD Hydraulic conductivity [L2/T]
#' @param c Hydraulic resistance [T]
#' @param h Constant in the definition of the stress (= sudden drawdown) [L]
#' @example calc_stress_response(x=100, t=1, S=0.001, kD=250, n=0, h=1)
#' @return A dataframe with the columns x, t, Phi, and h
calc_stress_response <- function(x, t, S, kD, c, h) {
      expand.grid(x = x, t = t) |>
            dplyr::rowwise() |>
            dplyr::mutate(
                  Labda = sqrt(c * kD),
                  Eta = 1 / (c * S),
                  h = h,
                  Phi = Br_123_32(x, t, h, Eta, Labda)
            ) |>
            dplyr::select(x, t, Phi, h)
}

#' Plot the Head change Phi [L] as a function of distance x [L] for different time points t
#'
#' @param df Dataframe with the results of calc_stress_response
#' @return A ggplot object
plot_stress_response_x_Phi <- function(df, title_text="Response surfacewater head changes") {
      p <- ggplot(df, aes(
            x = x,
            y = Phi,
            color = factor(t)
      )) +
            geom_line(linewidth = 1) +
            labs(
                  title = title_text,
                  x = "Distance [L]",
                  y = "Head change Phi [L]",
                  color = "Time [T]"
            ) +
            theme_minimal() +
            theme(
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16, margin = margin(t = 30)),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 16, margin = margin(b = 20)),
                  panel.border = element_rect(
                        color = "black",
                        fill = NA,
                        size = 0.5
                  ),
                  legend.box.background = element_rect(color = "black", size = 0.5)
            )
      # Convert ggplot to plotly
      plotly::ggplotly(p)
}

#' Plot the Head change Phi [L] as a function of time t [T] at different distances x [L]
#'
#' @param df Dataframe with the results of calc_stress_response
#' @return A ggplot object
plot_stress_response_t_Phi <- function(df, t_h = NULL, title_text="Response surfacewater head changes") {
      p <- ggplot(df, aes(
            x = t,
            y = Phi,
            color = factor(x)
      )) +
            geom_line(linewidth = 1) +
            labs(
                  title = title_text,
                  x = "Time [T]",
                  y = "Head change Phi [L]",
                  color = "Distance [L]"
            ) +
            theme_minimal() +
            theme(
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16, margin = margin(t = 30)),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 16, margin = margin(b = 20)),
                  panel.border = element_rect(
                        color = "black",
                        fill = NA,
                        size = 0.5
                  ),
                  legend.box.background = element_rect(color = "black", size = 0.5)
            )
      # Convert ggplot to plotly
      plotly::ggplotly(p)
}

############# Function related to simulation

#' Calculate the stress period index i, the stress a and the time t to the next output time for multiple values of t
#' 
#' @param t Vector of time values [T]
#' @param t_h Data frame with columns t [T] and h [L]
#' @example stress_period(t = c(1, 5, 10), t_h = data.frame(t = c(0, 2, 4), h = c(1, -1)))
#' @return A dataframe with columns t, sp, h, t0, and dt
stress_period <- function(t, t_h) {
      result_list <- lapply(t, function(single_t) {
            i <- which(t_h$t < single_t)
            data.frame(t = single_t, sp = i, h = t_h$h[i], t0 = t_h$t[i], dt = single_t - t_h$t[i])
      })
      result_df <- do.call(rbind, result_list)
      return(result_df)
}

#' Calculate the stress period data frame extended with a single value of x and multiple values of t
#' 
#' @param x Single value of distance [L]
#' @param t Vector of time values [T]
#' @param t_h Data frame with columns t [T] and h [L]
#' @example distance_stress_period_single(x = 5, t = c(1, 5, 10), t_h = data.frame(t = c(0, 2, 4), h = c(1, -1)))
#' @return A dataframe with columns x, t, sp, h, t0, and dt
distance_stress_period_single <- function(x, t, t_h) {
      data.frame(x = x, stress_period(t, t_h))
}

#' Calculate the stress period data frame extended with multiple values of x and multiple values of t
#' 
#' @param x Vector of distance values [L]
#' @param t Vector of time values [T]
#' @param t_h Data frame with columns t [T] and h [L]
#' @example distance_stress_period(x = c(5, 10, 15), t = c(1, 5, 10), t_h = data.frame(t = c(0, 2, 4), h = c(1, -1)))
#' @return A dataframe with columns x, t, sp, h, t0, and dt
distance_stress_period <- function(x, t, t_h) {
      result_list <- lapply(x, function(single_x) distance_stress_period_single(single_x, t, t_h))
      result_df <- do.call(rbind, result_list)
      return(result_df)
}

#' Simulate the stress response for given distances, output times, and stress-time series t_h
#' 
#' @param x Vector of distance values [L]
#' @param t Vector of output times [T]
#' @param S Storage coefficient [-]
#' @param kD Hydraulic conductivity [L2/T]
#' @param c Hydraulic resistance [T]
#' @param t_h Data frame with columns t [T] and h [L]
#' @example simulate_stress_response(x = c(5, 10, 100), t = c(1, 10, 100), S = 0.15, kD = 250, c=100, t_h = data.frame(t = c(0, 50), h = c(1, -1)))
#' @return A dataframe with summarized stress response
simulate_stress_response <- function(x, t, S, kD, c, t_h) {
      t_h <- t_h[order(t_h$t), ] # Order on time t
      t_h <- t_h[c(TRUE, diff(t_h$h) != 0), ] # Remove duplicates
      t_h$h[2:nrow(t_h)] <- diff(t_h$h) # Calculate differences
      df <- distance_stress_period(x, t, t_h)
      df_stress_response <- data.frame(b = mapply(calc_stress_response, x=df$x, t=df$dt, h=df$h, MoreArgs = list(S = S, kD = kD, c=c))) |>
            t() |> as.data.frame() |> dplyr::rename(dt = t)
      row.names(df_stress_response) <- NULL
      df_stress_response$t <- df$t
      df_stress_response <- data.frame(x = unlist(df_stress_response$x), t = unlist(df_stress_response$t), Phi = unlist(df_stress_response$Phi)) |>
            dplyr::group_by(x, t) |> dplyr::summarise(Phi = sum(Phi))
      return(df_stress_response)
}

# Define UI 
ui <- fluidPage(

    # Application title
    titlePanel("Stress Response Plot"),
    tabsetPanel(
          tabPanel("System definition",
                   sidebarLayout(
                         sidebarPanel(
                               numericInput("h", "Drawdown of surf. water level h [L]:", -1),
                               numericInput("kD", "Hydraulic conductivity kD [L2/T]:", value=2000, min=0.001),
                               numericInput("c", "Hydraulic resistance c [T]:", value=2000, min=10),
                               numericInput("S", "Storage coefficient S x 10-5 [-]:", value=0.5, min=0),
                               
                               bsTooltip("h", "value <> 0 [L]", "top", options = list(container = "body")),
                               bsTooltip("kD", "value > 0 [L2/T]", "top", options = list(container = "body")),
                               bsTooltip("c", "value > 10 [T]", "top", options = list(container = "body")),
                               bsTooltip("S x 10-5", "value > 0 [-]", "top", options = list(container = "body")),
                               
                               downloadButton("downloadData", "Download input data"), br(), br(),
                               
                               fileInput("uploadData", "Upload input data", accept = c(".csv")),
                               
                               tags$a(href = "https://github.com/KeesVanImmerzeel/Brug1DLeaky/tree/master", "Documentation")
                         ),
                         mainPanel(
                               
                         )
                   )
          ),
          tabPanel("System Plots", tabsetPanel(
                tabPanel("x, result",
                         sidebarLayout(
                               sidebarPanel(
                                     actionButton("plot_button_x", "Refresh", class = "btn-warning"),
                                     br(), br(),
                                     
                                     numericInput("num_points_x", "Number of points in x-array:", value = 100, min = 1),
                                     numericInput("min_x", "Minimum value of x:", value = 1, min = 0),
                                     numericInput("max_x", "Maximum value of x:", value = 1000, min = 0),
                                     numericInput("num_points_t", "Number of points in t-array:", value = 5, min = 1, max = 10),
                                     numericInput("min_t", "Minimum value of t:", value = 1, min = 0.001),
                                     numericInput("max_t", "Maximum value of t:", value = 1000, min = 0),
                                     
                                     bsTooltip("num_points_x", "value >= 1", "top", options = list(container = "body")),
                                     bsTooltip("num_points_t", "1 <= value <= 10", "top", options = list(container = "body")),
                                     bsTooltip("min_x", "value >= 0", "top", options = list(container = "body")),
                                     bsTooltip("max_x", "value > min_x", "top", options = list(container = "body")),
                                     bsTooltip("min_t", "value > 0", "top", options = list(container = "body")),
                                     bsTooltip("max_t", "value >= min_t", "top", options = list(container = "body")),
                                     
                                     br(),
                                     downloadButton("downloadResults", "Download Results")
                               ),
                               mainPanel(
                                     plotly::plotlyOutput("stressPlotS"), br(), br(),
                                     tableOutput("resultsTable")
                               )
                         )
                ),
                tabPanel("t, result",
                         sidebarLayout(
                               sidebarPanel(
                                     actionButton("plot_button_t", "Refresh", class = "btn-warning"),
                                     br(), br(),
                                     
                                     numericInput("num_points_x_t", "Number of points in x-array:", value = 6, min = 1, max=10),
                                     numericInput("min_x_t", "Minimum value of x:", value = 0, min = 0),
                                     numericInput("max_x_t", "Maximum value of x:", value = 1000, min = 0),
                                     numericInput("num_points_t_t", "Number of points in t-array:", value = 100, min = 1),
                                     numericInput("min_t_t", "Minimum value of t:", value = 1, min = 0.001),
                                     numericInput("max_t_t", "Maximum value of t:", value = 1000, min = 0),
                                     
                                     bsTooltip("num_points_x_t","1 <= value <= 10" , "top", options = list(container = "body")),
                                     bsTooltip("num_points_t_t", "value >= 1", "top", options = list(container = "body")),
                                     bsTooltip("min_x_t", "value > 0", "top", options = list(container = "body")),
                                     bsTooltip("max_x_t", "value > min_x_t", "top", options = list(container = "body")),
                                     bsTooltip("min_t_t", "value >= 0", "top", options = list(container = "body")),
                                     bsTooltip("max_t_t", "value >= min_t_t", "top", options = list(container = "body")),
                                     
                                     br(),
                                     downloadButton("downloadResults_t", "Download Results")
                               ),
                               mainPanel(
                                     plotly::plotlyOutput("stressPlotS_t"), br(), br(),
                                     tableOutput("resultsTable_t")
                               )
                         )
                )
          )),
          tabPanel("Simulation definition",
                   sidebarLayout(
                         sidebarPanel(
                               fileInput("file1", "Upload Spreadsheet",
                                         accept = c(".csv", ".xlsx"))
                         ),
                         mainPanel(
                               plotly::plotlyOutput("stress_plot"),
                               tableOutput("stress_sequence")
                         )
                   )
          ), 
          tabPanel("Simulation plots", tabsetPanel(
                tabPanel("x, result",
                         sidebarLayout(
                               sidebarPanel(
                                     actionButton("sim_plot_button_x", "Refresh", class = "btn-warning"),
                                     br(), br(),
                                     
                                     numericInput("sim_num_points_x", "Number of points in x-array:", value = 100, min = 1),
                                     numericInput("sim_min_x", "Minimum value of x:", value = 1, min = 0),
                                     numericInput("sim_max_x", "Maximum value of x:", value = 1000, min = 0),
                                     numericInput("sim_num_points_t", "Number of points in t-array:", value = 7, min = 1, max = 10),
                                     numericInput("sim_min_t", "Minimum value of t:", value = 7, min = 0.001),
                                     numericInput("sim_max_t", "Maximum value of t:", value = 120, min = 0),
                                     
                                     bsTooltip("sim_num_points_x", "value >= 1", "top", options = list(container = "body")),
                                     bsTooltip("sim_num_points_t", "1 <= value <= 10", "top", options = list(container = "body")),
                                     bsTooltip("sim_min_x", "value >= 0", "top", options = list(container = "body")),
                                     bsTooltip("sim_max_x", "value > min_x", "top", options = list(container = "body")),
                                     bsTooltip("sim_min_t", "value > 0", "top", options = list(container = "body")),
                                     bsTooltip("sim_max_t", "value >= min_t", "top", options = list(container = "body")),
                                     
                                     br(),
                                     downloadButton("sim_downloadResults", "Download Results")
                               ),
                               mainPanel(
                                     plotly::plotlyOutput("sim_stressPlotS"), br(), br(),
                                     tableOutput("sim_resultsTable")
                               )
                         )
                ),
                tabPanel("t, result",
                         sidebarLayout(
                               sidebarPanel(
                                     actionButton("sim_plot_button_t", "Refresh", class = "btn-warning"),
                                     br(), br(),
                                     
                                     numericInput("sim_num_points_x_t", "Number of points in x-array:", value = 6, min = 1, max=10),
                                     numericInput("sim_min_x_t", "Minimum value of x:", value = 100, min = 0),
                                     numericInput("sim_max_x_t", "Maximum value of x:", value = 1000, min = 0),
                                     numericInput("sim_num_points_t_t", "Number of points in t-array:", value = 100, min = 1),
                                     numericInput("sim_min_t_t", "Minimum value of t:", value = 10, min = 0.001),
                                     numericInput("sim_max_t_t", "Maximum value of t:", value = 150, min = 0),
                                     
                                     bsTooltip("sim_num_points_x_t","1 <= value <= 10" , "top", options = list(container = "body")),
                                     bsTooltip("sim_num_points_t_t", "value >= 1", "top", options = list(container = "body")),
                                     bsTooltip("sim_min_x_t", "value > 0", "top", options = list(container = "body")),
                                     bsTooltip("sim_max_x_t", "value > min_x_t", "top", options = list(container = "body")),
                                     bsTooltip("sim_min_t_t", "value >= 0", "top", options = list(container = "body")),
                                     bsTooltip("sim_max_t_t", "value >= min_t_t", "top", options = list(container = "body")),
                                     
                                     br(),
                                     downloadButton("sim_downloadResults_t", "Download Results")
                               ),
                               mainPanel(
                                     plotly::plotlyOutput("sim_stressPlotS_t"), br(), br(),
                                     tableOutput("sim_resultsTable_t")
                               )
                         )
                )
          ))
    )
)

# Server logic.
server <- function(input, output, session) {
      result_x <- reactiveVal()
      result_t <- reactiveVal()
      sim_result_x <- reactiveVal()
      sim_result_t <- reactiveVal()
      
      observeEvent(input$plot_button_x, {
            x_values <- seq(input$min_x, input$max_x, length.out = input$num_points_x)
            t_values <- round(pracma::logseq(input$min_t, input$max_t, input$num_points_t),1)
            result_x(calc_stress_response(x = x_values, t = t_values, S = input$S*10^-5, kD = input$kD, c = input$c, h = input$h))
            output$stressPlotS <- renderPlotly({
                  plot_stress_response_x_Phi(result_x(), title_text= paste("Sudden drawdown of the surface water level with h=", input$h, "[L]")) 
            })
            output$resultsTable <- renderTable({
                  result_x()
            })
      })
      
      observeEvent(input$plot_button_t, {
            x_values <- seq(input$min_x_t, input$max_x_t, length.out = input$num_points_x_t)
            t_values <- seq(input$min_t_t, input$max_t_t, length.out = input$num_points_t_t)
            
            result_t(calc_stress_response(x = x_values, t = t_values, S = input$S*10^-5, kD = input$kD, c = input$c, h = input$h))
            title_text <- paste("Sudden drawdown of the surface water level with h=", input$h, "[L]")
            output$stressPlotS_t <- renderPlotly({
                  plot_stress_response_t_Phi(result_t(), title_text= paste("Sudden drawdown of the surface water level with h=", input$h, "[L]"))
            })
            output$resultsTable_t <- renderTable({
                  result_t()
            })
      })  
      
      
      output$downloadData <- downloadHandler(
            filename = function() {
                  paste("Brug1DLeaky_input_data_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  input_data <- data.frame(
                        h = input$h,
                        kD = input$kD,
                        c = input$c,                        
                        S = input$S,
                        
                        num_points_x = input$num_points_x,
                        min_x = input$min_x,
                        max_x = input$max_x,
                        num_points_t = input$num_points_t,
                        min_t = input$min_t,
                        max_t = input$max_t,
                        
                        num_points_x_t = input$num_points_x_t,
                        min_x_t = input$min_x_t,
                        max_x_t = input$max_x_t,
                        num_points_t_t = input$num_points_t_t,
                        min_t_t = input$min_t_t,
                        max_t_t = input$max_t_t,
                        
                        sim_num_points_x = input$sim_num_points_x,
                        sim_min_x = input$sim_min_x,
                        sim_max_x = input$sim_max_x,
                        sim_num_points_t = input$sim_num_points_t,
                        sim_min_t = input$sim_min_t,
                        sim_max_t = input$sim_max_t,
                        
                        sim_num_points_x_t = input$sim_num_points_x_t,
                        sim_min_x_t = input$sim_min_x_t,
                        sim_max_x_t = input$sim_max_x_t,
                        sim_num_points_t_t = input$sim_num_points_t_t,
                        sim_min_t_t = input$sim_min_t_t,
                        sim_max_t_t = input$sim_max_t_t
                  )
                  write.csv2(input_data, file, quote = FALSE, row.names = FALSE)
            }
      )
      
      output$downloadResults <- downloadHandler(
            filename = function() {
                  paste("Brug1DLeaky_x_vs_results_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  write.csv2(result_x(), file, quote = FALSE, row.names = FALSE)
            }
      )
      
      output$downloadResults_t <- downloadHandler(
            filename = function() {
                  paste("Brug1DLeaky_t_vs_results_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  write.csv2(result_t(), file, quote = FALSE, row.names = FALSE)
            }
      )
      
      output$sim_downloadResults <- downloadHandler(
            filename = function() {
                  paste("Brug1DLeaky_x_vs_sim_results_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  write.csv2(sim_result_x(), file, quote = FALSE, row.names = FALSE)
            }
      )
      
      output$sim_downloadResults_t <- downloadHandler(
            filename = function() {
                  paste("Brug1DLeaky_t_vs_sim_results_", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                  write.csv2(sim_result_t(), file, quote = FALSE, row.names = FALSE)
            }
      )
      
      observeEvent(input$uploadData, {
            req(input$uploadData)
            input_data <- read.csv(input$uploadData$datapath)
            updateNumericInput(session, "h", value = input_data$h)    
            updateNumericInput(session, "kD", value = input_data$kD)            
            updateRadioButtons(session, "c", selected = input_data$c)
            updateNumericInput(session, "S", value = input_data$S)
            
            updateNumericInput(session, "num_points_x", value = input_data$num_points_x)
            updateNumericInput(session, "min_x", value = input_data$min_x)
            updateNumericInput(session, "max_x", value = input_data$max_x)
            updateNumericInput(session, "num_points_t", value = input_data$num_points_t)
            updateNumericInput(session, "min_t", value = input_data$min_t)
            updateNumericInput(session, "max_t", value = input_data$max_t)
            
            updateNumericInput(session, "num_points_x_t", value = input_data$num_points_x_t)
            updateNumericInput(session, "min_x_t", value = input_data$min_x_t)
            updateNumericInput(session, "max_x_t", value = input_data$max_x_t)
            updateNumericInput(session, "num_points_t_t", value = input_data$num_points_t_t)
            updateNumericInput(session, "min_t_t", value = input_data$min_t_t)
            updateNumericInput(session, "max_t_t", value = input_data$max_t_t)
            
            updateNumericInput(session, "sim_num_points_x", value = input_data$sim_num_points_x)
            updateNumericInput(session, "sim_min_x", value = input_data$sim_min_x)
            updateNumericInput(session, "sim_max_x", value = input_data$sim_max_x)
            updateNumericInput(session, "sim_num_points_t", value = input_data$sim_num_points_t)
            updateNumericInput(session, "sim_min_t", value = input_data$sim_min_t)
            updateNumericInput(session, "sim_max_t", value = input_data$sim_max_t)
            
            updateNumericInput(session, "sim_num_points_x_t", value = input_data$sim_num_points_x_t)
            updateNumericInput(session, "sim_min_x_t", value = input_data$sim_min_x_t)
            updateNumericInput(session, "sim_max_x_t", value = input_data$sim_max_x_t)
            updateNumericInput(session, "sim_num_points_t_t", value = input_data$sim_num_points_t_t)
            updateNumericInput(session, "sim_min_t_t", value = input_data$sim_min_t_t)
            updateNumericInput(session, "sim_max_t_t", value = input_data$sim_max_t_t)
      })
      
      
      ######### Simulate in- and output
      
      stress_sequence <- reactive({
            req(input$file1)
            inFile <- input$file1
            
            if (grepl("\\.csv$", inFile$name)) {
                  df <- read.csv(inFile$datapath)
            } else if (grepl("\\.xlsx$", inFile$name)) {
                  df <- readxl::read_excel(inFile$datapath)
            }
            
            #df <- df[, c("Time", "h")]
            names(df) <- c("Time", "h")
            df <- df[order(df$Time), ]
            #df <- df[c(TRUE, diff(df$h) != 0), ] # Remove duplicates
            return(df)
      })
            
      output$stress_plot <- renderPlotly({
            req(stress_sequence())
            df <- stress_sequence()
            
            # Ensure values remain the same until the next timestep
            df <- df %>%
                  mutate(Time_end = lead(Time, default = max(Time)),
                         h_end = h) %>%
                  tidyr::pivot_longer(cols = c("Time", "Time_end"), names_to = "Time_type", values_to = "Time_value") %>%
                  arrange(Time_value)
            
            ylab <- "h [L]"
            p <- ggplot(df, aes(x = Time_value, y = h)) +
                  geom_step(linewidth = 1) +
                  labs(title = "Stress Sequence", x = "Time [T]", y = ylab) +
                  theme(
                        axis.text = element_text(size = 14),
                        axis.title = element_text(size = 16, margin = margin(t = 30)),
                        legend.text = element_text(size = 14),
                        legend.title = element_text(size = 16, margin = margin(b = 20)),
                        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                        legend.box.background = element_rect(color = "black", size = 0.5)
                  )
            # Convert ggplot to plotly
            plotly::ggplotly(p)
            
      })      
      
      output$stress_sequence <- renderTable({
            req(stress_sequence())
            stress_sequence()
      })      
            
      observeEvent(input$sim_plot_button_x, {
            x_values <- seq(input$sim_min_x, input$sim_max_x, length.out = input$sim_num_points_x)
            t_values <- round(pracma::logseq(input$sim_min_t, input$sim_max_t, input$sim_num_points_t),1)
            
            t_h <- stress_sequence()
            names(t_h) <- c("t", "h")
            sim_result_x(simulate_stress_response(x = x_values, t = t_values, S = input$S*10^-5, kD = input$kD, c = input$c, t_h = t_h ))
            output$sim_stressPlotS <- renderPlotly({
                  plot_stress_response_x_Phi(sim_result_x())
            })
            output$sim_resultsTable <- renderTable({
                  sim_result_x()
            })
      })

      observeEvent(input$sim_plot_button_t, {
            x_values <- seq(input$sim_min_x_t, input$sim_max_x_t, length.out = input$sim_num_points_x_t)
            t_values <- seq(input$sim_min_t_t, input$sim_max_t_t, length.out = input$sim_num_points_t_t)
            
            t_h <- stress_sequence()
            names(t_h) <- c("t", "h")
            sim_result_t(simulate_stress_response(x = x_values, t = t_values, S = input$S*10^-5, kD = input$kD, c = input$c, t_h = t_h))
            output$sim_stressPlotS_t <- renderPlotly({
                  plot_stress_response_t_Phi(sim_result_t(), t_a)
            })
            output$sim_resultsTable_t <- renderTable({
                  sim_result_t()
            })
      })            
}

# Run the application 
shinyApp(ui = ui, server = server)
