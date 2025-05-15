# Load libraries yang diperlukan
library(shiny)
library(forecast)
library(ggplot2)
library(readr)
library(tseries)
library(dplyr)
library(bslib)
library(DT)

#Tampilan UI
ui <- fluidPage(
  # Tema Bootstrap 5 dengan font Google Poppins
  theme = bs_theme(
    version = 5,
    bootswatch = "sketchy", # Tema awal 
    base_font = font_google("Poppins"),
    heading_font = font_google("Poppins"),
    primary = "#007BFF"
  ),
  
  # Judul aplikasi
  titlePanel("📊 Forecasting Triple Exponential Smoothing (Holt-Winters)"),
  
  # Layout halaman dibagi menjadi sidebar dan main panel
  sidebarLayout(
    sidebarPanel(
      # Pilihan sumber data
      radioButtons("data_choice", "Pilih Sumber Data:",
                   choices = c("Unggah Data" = "upload", 
                               "AirPassengers" = "airpassengers",
                               "UKgas" = "ukgas",
                               "Nottem (Suhu)" = "nottem")),
      
      # Input file CSV hanya ditampilkan jika pilihan adalah 'Unggah Data'
      conditionalPanel(
        condition = "input.data_choice == 'upload'",
        fileInput("file", "Unggah File CSV", accept = ".csv"),
        helpText("Format: kolom 1 = Tanggal (YYYY-MM-DD), kolom 2 = Nilai")
      ),
      
      # Pilihan jenis musiman
      selectInput("seasonal_type", "Tipe Musiman:", 
                  choices = c("Multiplicative" = "multiplicative", 
                              "Additive" = "additive")),
      
      # Input periode musiman dan periode prediksi
      numericInput("seasonal_period", "Periode Musiman:", value = 0, min = 1),
      numericInput("horizon", "Periode Prediksi:", value = 0, min = 1),
      
      # Tombol untuk memulai analisis
      actionButton("analyze", "🔍 Analisis & Forecast"),
      br(), br(),
      
      # Tombol untuk mengunduh hasil forecasting
      downloadButton("download_forecast", "⬇️ Unduh Forecast (.csv)")
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        # Tab untuk preview data
        tabPanel("📄 Preview Data",
                 h4("📥 Data Input"),
                 dataTableOutput("uploaded_data")
        ),
        
        # Tab untuk analisis data
        tabPanel("📈 Analisis Data",
                 h4("📊 Hasil Analisis"),
                 verbatimTextOutput("summary"),         # Statistik ringkasan
                 plotOutput("plot"),                    # Plot data asli
                 plotOutput("seasonal_plot"),           # Plot dekomposisi musiman
                 verbatimTextOutput("seasonality_check"),  # Cek pola musiman
                 verbatimTextOutput("stationarity_check")  # Cek stasioneritas
        ),
        
        # Tab untuk forecasting
        tabPanel("🔮 Forecasting",
                 h4("📈 Hasil Forecasting"),
                 verbatimTextOutput("model_params"),    # Parameter model Holt-Winters
                 verbatimTextOutput("forecast_result"), # Output forecast lengkap
                 verbatimTextOutput("mape_result"),     # MAPE dari hasil fitting
                 plotOutput("forecast_plot")            # Plot forecasting
        )
      )
    )
  )
)

# Server 
server <- function(input, output, session) {
  values <- reactiveValues(data = NULL, ts_data = NULL, forecast_tbl = NULL)
  observeEvent(input$analyze, {
    
    # Bagian pemrosesan data
    if (input$data_choice == "upload") {
      # Validasi file yang diunggah
      req(input$file)
      raw_data <- read_csv2(input$file$datapath, show_col_types = FALSE)
      data <- raw_data[, 1:2]
      colnames(data) <- c("Tanggal", "Nilai")
      data$Tanggal <- as.Date(data$Tanggal, format = "%Y-%m-%d")
      data <- data %>% arrange(Tanggal)
      
      # Konversi ke time series
      start_year <- as.numeric(format(min(data$Tanggal), "%Y"))
      start_month <- as.numeric(format(min(data$Tanggal), "%m"))
      ts_data <- ts(data$Nilai, frequency = input$seasonal_period, start = c(start_year, start_month))
    } else {
      # Ambil data default dari R
      ts_data <- switch(input$data_choice,
                        "airpassengers" = AirPassengers,
                        "ukgas" = UKgas,
                        "nottem" = nottem)
      data <- data.frame(Tanggal = time(ts_data), Nilai = as.numeric(ts_data))
    }
    
    values$data <- data
    values$ts_data <- ts_data
    
    # Tampilkan data ke UI
    output$uploaded_data <- renderDT({
      datatable(values$data,
                options = list(
                  columnDefs = list(
                    list(className = 'dt-center', targets = "_all"),
                    list(width = '100px', targets = "_all")
                  )
                ),
                rownames = FALSE,
                class = 'cell-border stripe hover compact'
      )
    })
    
    # Ringkasan statistik nilai
    output$summary <- renderPrint({
      summary(data$Nilai)
    })
    
    # Plot time series data asli
    output$plot <- renderPlot({
      ggplot(data, aes(x = Tanggal, y = Nilai)) +
        geom_line(color = "blue", linewidth = 1) +
        labs(title = "Grafik Time Series", x = "Tanggal", y = "Nilai") +
        theme_minimal()
    })
    
    # Dekomposisi musiman 
    decomposed <- stl(ts_data, s.window = "periodic")
    
    output$seasonal_plot <- renderPlot({
      autoplot(decomposed) +
        labs(title = "Dekomposisi Komponen") +
        theme_minimal()
    })
    
    # Cek apakah ada komponen musiman
    seasonal_present <- sum(decomposed$time.series[, "seasonal"]^2) > 0
    
    output$seasonality_check <- renderText({
      if (seasonal_present) {
        "✅ Data memiliki pola musiman. Lanjut ke uji stasioneritas."
      } else {
        "⚠️ Data tidak memiliki pola musiman."
      }
    })
    
    # Jika musiman, lanjut ke uji stasioneritas dan forecasting
    if (seasonal_present) {
      # Uji ADF
      adf_result <- adf.test(ts_data)
      output$stationarity_check <- renderText({
        pval <- round(adf_result$p.value, 4)
        if (pval > 0.05) {
          paste0("Data tidak stasioner (p-value = ", pval, " > 0.05).")
        } else {
          paste0("Data stasioner (p-value = ", pval, " <= 0.05).")
        }
      })
      
      # Model Holt-Winters
      model <- HoltWinters(ts_data, seasonal = input$seasonal_type)
      forecast_result <- forecast(model, h = input$horizon)
      
      # Parameter model
      output$model_params <- renderPrint({
        list(
          alpha = model$alpha,
          beta = model$beta,
          gamma = model$gamma,
          catatan = ifelse(adf_result$p.value > 0.05,
                           "✅ Cocok dengan Holt-Winters karena data tidak stasioner.",
                           "⚠️ Holt-Winters tetap bisa digunakan meskipun data stasioner.")
        )
      })
      
      # Output forecasting
      output$forecast_result <- renderPrint({ forecast_result })
      
      # Hitung MAPE
      fitted_vals <- fitted(model)[,1]
      original_vals <- ts_data[1:length(fitted_vals)]
      mape <- mean(abs((original_vals - fitted_vals) / original_vals)) * 100
      
      output$mape_result <- renderPrint({
        paste0("MAPE: ", round(mape, 2), "%")
      })
      
      # Simpan hasil forecasting sebagai data.frame
      forecast_df <- data.frame(
        Periode = time(forecast_result$mean),
        Forecast = as.numeric(forecast_result$mean),
        Lower80 = as.numeric(forecast_result$lower[,1]),
        Upper80 = as.numeric(forecast_result$upper[,1])
      )
      values$forecast_tbl <- forecast_df
      
      # Plot hasil forecasting
      output$forecast_plot <- renderPlot({
        autoplot(forecast_result) +
          labs(title = "Hasil Forecast", x = "Periode", y = "Nilai") +
          theme_minimal()
      })
    }
  })
  
  # Fungsi unduh hasil forecast
  output$download_forecast <- downloadHandler(
    filename = function() {
      paste0("forecast_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$forecast_tbl)
      write.csv(values$forecast_tbl, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
