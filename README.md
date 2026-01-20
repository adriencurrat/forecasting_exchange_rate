# Forecasting the Swiss nominal effective exchange rate

**Authors:** Isaac Graber, Mathéo Bourgeois, Adrien Currat  
**Institution:** University of Lausanne, Msc. in Economics  
**Course:** Economic Forecasting for Decision-Making  

## Introduction
This project replicates a real-world forecasting exercise conducted in the context of the Swiss National Bank (SNB). Acting as a macroeconomic forecasting team, we forecasted the evolution of the Swiss nominal effective exchange rate (NEER) of Switzerland up to 2026.

Using a theoretical framework based on the Taylor rule and the uncovered interest parity (UIP) condition, we assessed the impact of inflation and output gap differentials on bilateral exchange rates with major trading partners (US, EU, China, and UK), weighted by trade importance. The results showed that the model delivered accurate forecasts and outperformed a random walk.

This project received the maximum grade in the course. The work was also highlighted by Professor Grobéty in a public [LinkedIn post](https://www.linkedin.com/feed/update/urn:li:activity:7275406513182060545/), acknowledging the quality of the analysis.

## Project Structure

### 1. Project workflow
To run the project:

1. Place all scripts and datasets in the same working directory.  
2. Open the `Main.R` script in RStudio.  
3. Ensure the working directory is correctly set in `Main.R`.  
4. Source the `Main.R` script to execute all steps in the correct order.

### 2. Scripts overview

#### 3.1 Packages Loading.R
**Purpose:** Ensures that all required R packages are installed and loaded. The script sets the CRAN mirror and installs missing packages automatically.

#### 3.2 Data Treatment.R
**Purpose:** Handles data cleaning, preprocessing, and interpolation. Raw macroeconomic data are transformed into usable formats for the analysis.

#### 3.3 Data Analysis.R
**Purpose:** Conducts stationarity testing, seasonal adjustment, and exploratory data analysis. Outputs include summary statistics and graphical representations of key variables.

#### 3.4 Model Estimation.R
**Purpose:** Estimates Taylor rule and Uncovered Interest Parity (UIP)–based models for exchange rate forecasting. Outputs include estimated coefficients and diagnostic metrics.

#### 3.5 Point Forecasting.R
**Purpose:** Generates point forecasts for NEER changes at multiple horizons (1, 6, 12, and 18 months).

#### 3.6 Forecast Evaluation.R
**Purpose:** Evaluates forecast accuracy using standard metrics and compares model performance against a random walk benchmark.

#### 3.7 Density Forecast.R
**Purpose:** Assesses the quality of density forecasts and prediction intervals. Includes Likelihood Ratio (LR) tests, Probability Integral Transform (PIT) analysis, and autocorrelation tests.

