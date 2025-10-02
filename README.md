# Option Pricing and Hedging under Additive Models

## Overview
The project focuses on calibrating **Additive Bachelier** and **Additive Logistic** models to market option data during the extraordinary conditions of April 2020, when oil prices briefly turned negative.  
The analysis is extended to **exotic option pricing** and **hedging strategies**, with implementations in both **MATLAB** and **Python**.

---

## Features
- **Discount Factors and Market Stress**
  - Recovery of discount factors and forward prices directly from option data.  
  - Analysis of April 2020 market dislocations.  

- **Additive Bachelier Model**
  - Calibration of ATM volatilities.  
  - Estimation of skewness (η) and vol-of-vol (κ) parameters.  
  - Volatility surface fitting via Fourier (Lewis formula).  

- **Additive Logistic Model**
  - Calibration of σ and H parameters.  
  - Comparison with empirical ATM volatility curves.  
  - Self-similar additive process for option pricing.  

- **Exotic Option Pricing**
  - Monte Carlo pricing of a binary (digital) option with average-forward payoff.  
  - Comparison of results under both models.  

- **Hedging Strategies**
  - Vega hedging via ATM volatility bumping & Pareto selection.  
  - Hedging sensitivities to η (bull spread) and κ (strangle).  
  - Delta hedging with futures.  
  - P&L analysis of hedged vs unhedged portfolios.  

- **Dual Implementation**
  - **MATLAB**: main calibration and hedging routines.  
  - **Python**: translation of core modules, validation of consistency.  

---

## Project Structure
```
├── matlab/                         # MATLAB implementation
│   ├── calibration_bachelier.m
│   ├── calibration_logistic.m
│   ├── exotic_pricing.m
│   ├── hedging_strategies.m
│   └── utils/
│
├── python/                         # Python implementation
│   ├── Calibration_Bachelier.ipynb
│   ├── Calibration_Logistic.ipynb
│   ├── Exotic_Option_Pricing.ipynb
│   └── Hedging_Strategies.ipynb
│
├── data/                           # Market option data (WTI options, June 2020)
│   ├── datacalls/
│   └── dataputs/
│
├── report/                         # Final report (PDF)
│   └── Assignment2_Report.pdf
│
└── README.md                       # Project documentation
```

---

## Results
- The **Additive Bachelier model** achieved excellent fit for short maturities and captured left skew in implied volatilities.  
- The **Additive Logistic model** provided a flexible alternative, closely matching ATM volatility term structure.  
- **Exotic option pricing** yielded consistent results across both models (difference < 0.03%).  
- **Hedging** improved portfolio stability significantly:  
  - Unhedged P&L ≈ –220k USD after 1 week.  
  - Hedged P&L ≈ +180k USD after 1 week.  

---

## Technologies
- **MATLAB R2023b** – main implementation.  
- **Python 3.x** – validation and alternative simulations.  
- Libraries: `numpy`, `scipy`, `matplotlib`, `sympy` for characteristic functions.  
