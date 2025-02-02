Here's a **README.md** file incorporating the purpose of your repository, a brief explanation of the scripts, installation instructions, usage, and the included figures.

---

## 🔥 Hydrogen Flame Speed Model

This repository contains Python scripts for simulating **laminar and turbulent flame speeds** in **hydrogen-air combustion** using **Cantera**. The models compare different mechanisms for calculating laminar flame speed and use a **flame stretch model** for predicting turbulent flame speed.

---

### 📌 Features
- Simulates **laminar flame speed** using different reaction mechanisms (H2O2 and GRI 3.0).
- Implements a **flame stretch model** to compute **turbulent flame speed**.
- Provides **experimental data comparison** for validation.
- Saves results as CSV files and generates **high-resolution plots**.
- Uses **Cantera** for chemical kinetics calculations.

---

## 📂 Repository Structure
```
📁 Hydrogen_Flame_Speed_Model
│── LFS.py                  # Laminar Flame Speed Simulation
│── TFS_EDC_Model.py         # Turbulent Flame Speed Model using EDC
│── LFS_reference.csv        # Experimental Data for Laminar Flame Speed
│── README.md                # Documentation
│── flame_speed_comparison_4k.png  # Plot of flame speeds comparison
│── mechanism_comparison.png       # Mechanism comparison results
```

---

## 🛠️ Installation

Ensure you have **Python 3.8+** installed. Then, install the required dependencies:

```sh
pip install numpy pandas scipy matplotlib cantera
```

---

## 🚀 Usage

### **1️⃣ Running the Laminar Flame Speed Model**
To simulate laminar flame speed:

```sh
python LFS.py
```

This script:
- Computes **laminar flame speed** for different equivalence ratios.
- Compares results from **H2O2 and GRI 3.0 mechanisms**.
- Saves results in CSV format.
- Generates a plot for visualization.

### **2️⃣ Running the Turbulent Flame Speed Model**
To compute turbulent flame speed:

```sh
python TFS_EDC_Model.py
```

This script:
- Uses the computed **laminar flame speed**.
- Applies a **flame stretch correction model**.
- Compares computed results against **experimental data**.
- Saves results as a CSV file.

---

## 📊 Results

### **Laminar & Turbulent Flame Speed Comparison**
![Flame Speed Comparison](flame_speed_comparison_4k.png)

- The **blue curve** represents the **laminar flame speed**.
- The **red curve** represents the **computed turbulent flame speed**.
- The **black stars** represent **experimental data**.

---

### **Mechanism Comparison for Laminar Flame Speed**
![Mechanism Comparison](mechanism_comparison.png)

- The **black line** represents the **H2O2 mechanism**.
- The **red dashed line** represents the **GRI 3.0 mechanism**.
- The **black circles** are **experimental data**.

---

## 📜 References
- **Magnussen, B. F. (1982)**: *A Simple Eddy Dissipation Concept Model for Turbulent Combustion*.
- **Burke, E. M., Singlitico, A., Morones, A. E., Guthe, F., Speth, L., & Monaghan, R. F. D. (2015)**. Progress towards a validated Cantera-based turbulent flame speed solver. Seventh Eur. Combust. Meet, 1-6.
- **Cantera Documentation**: https://cantera.org/

---

## 🏷 License
This project is released under the **MIT License**.
