# master-thesis-longitudinal-analysis
Longitudinal Analysis of Metabolic Parameters of Statin Therapy  in Diabetic Patients using GEE, LME, and Deep Learning models.

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0%2B-276DC3.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

本專案為碩士論文研究，基於真實世界電子健康紀錄（EHR），探討糖尿病病患服用 Statin 藥物與代謝指標之縱向關聯，並系統性比較傳統統計模型與時間序列深度學習模型的預測效能。

---

## 📊 專案亮點與資料規模
* **臨床真實數據**：納入 **1,790 位病患與 17,252 筆縱向觀測** 的大型追蹤資料。
* **雙軌建模方法**：結合傳統統計學（推論個體異質性）與現代深度學習（條件式預測）。
* **多指標評估**：針對 HbA1c、Glu-AC、LDL 與 TG 等關鍵代謝指標進行效能與可解釋性比較。

---

## 🧠 核心方法與架構
1. **資料前處理與特徵工程**（Python / R）：執行端到端（End-to-End）的資料清理、遺漏值處理與縱向資料塑形。
2. **傳統縱向統計模型**：
   * **GEE（廣義估計方程式）**
   * **LME（線性混合效應模型，含隨機截距與斜率）**
3. **時間序列深度學習模型**：
   * `LSTM-NARX`
   * `TCN` (Temporal Convolutional Network)
   * `ODE-LSTM` (簡化式常微分方程式 LSTM)
4. **模擬研究（Simulation Study）**：評估不同資料變異結構對模型預測誤差的影響。

---
